import os
import pathlib
import pkgutil
import shutil
import struct
import sys
import platform
import setuptools
import distutils.sysconfig as sysconfig
from distutils.core import setup
import distutils.command.build
from distutils.command.install_data import install_data
from subprocess import CalledProcessError, check_output, check_call

from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install_lib import install_lib
from setuptools.command.install_scripts import install_scripts
from setuptools.config.pyprojecttoml import read_configuration

try:
    from find_libpython import find_libpython

    p_lib = find_libpython()
    if os.environ.get('FLIP_SLASHES', False):
        p_lib = p_lib.replace('\\', '/')
    if p_lib is not None:
        python_lib_path = f'-DPYTHON_LIBRARY_PATH={p_lib}'
    else:
        python_lib_path = ''
except ImportError:
    p_dir = sysconfig.get_config_var("LIBDIR")
    p_lib = sysconfig.get_config_var("INSTSONAME")
    if p_dir is None or p_lib is None:
        python_lib_path = ''
    else:
        python_lib_path = f'-DPYTHON_LIBRARY_PATH={os.path.join(p_dir, p_lib)}'

try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel


    # this overrides standard naming of the wheel to not include
    # architecture or python dot version number

    class Bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            self.root_is_pure = False

        def get_tag(self):
            python, abi, plat = _bdist_wheel.get_tag(self)
            if platform.system() != 'Windows':
                python, abi = 'py3', 'none'
            return python, abi, plat
except ImportError:
    Bdist_wheel = None

BITS = struct.calcsize("P") * 8
PACKAGE_NAME = "crysfml_api"
SOURCE_DIR = '.'
COMPILER = 'gfortran'
if os.environ.get('FC', False):
    COMPILER = os.environ.get('FC')
print(f'Compiler set to: {COMPILER}')
CFML_OVERRIDE = []
if os.environ.get('CMAKE_GENERATOR', False):
    CFML_OVERRIDE += ["-G", f"{os.environ.get('CMAKE_GENERATOR')}"]

# We can use cmake provided from pip which (normally) gets installed at /bin
# Except that in the manylinux builds it's placed at /opt/python/[version]/bin/
# (as a symlink at least) which is *not* on the path.
# If cmake is a known module, import it and use it tell us its binary directory
if pkgutil.find_loader('cmake') is not None:
    import cmake

    CMAKE_BIN = cmake.CMAKE_BIN_DIR + os.path.sep + 'cmake'
else:
    CMAKE_BIN = 'cmake'


def get_cmake():
    return CMAKE_BIN


class CMakeExtension(Extension):
    """
    An extension to run the cmake build
    This simply overrides the base extension class so that setuptools
    doesn't try to build your sources for you
    """

    def __init__(self, name, sources=[]):
        super().__init__(name=name, sources=sources)


class InstallCMakeLibsData(install_data):
    """
    Just a wrapper to get the install data into the egg-info
    Listing the installed files in the egg-info guarantees that
    all of the package files will be uninstalled when the user
    uninstalls your package through pip
    """

    def run(self):
        """
        Out files are the libraries that were built using cmake
        """
        # There seems to be no other way to do this; I tried listing the
        # libraries during the execution of the InstallCMakeLibs.run() but
        # setuptools never tracked them, seems like setuptools wants to
        # track the libraries through package data more than anything...

        self.outfiles = self.distribution.data_files


class InstallCMakeLibs(install_lib):
    """
    Get the libraries from the parent distribution, use those as the outfiles
    Skip building anything; everything is already built, forward libraries to
    the installation step
    """

    def run(self):
        """
        Copy libraries from the bin directory and place them as appropriate
        """

        self.announce("Moving library files", level=3)
        # We have already built the libraries in the previous build_ext step
        self.skip_build = True

        # Depending on the files that are generated from your cmake
        # build chain, you may need to change the below code, such that
        # your files are moved to the appropriate location when the installation
        # is run

        # bin_dir = self.distribution.bin_dir
        bin_dir = os.path.abspath(os.path.join(self.distribution.bin_dir, '..'))

        libs = [os.path.join(bin_dir, _dir) for _dir in
                os.listdir(bin_dir) if
                os.path.isdir(os.path.join(bin_dir, _dir))]

        for lib in libs:
            shutil.move(lib, os.path.join(self.build_dir,
                                          os.path.basename(lib)))

        # Move the lib to the correct location.
        bin_dir = self.build_dir
        pyd_path = [os.path.join(bin_dir, _pyd) for _pyd in
                    os.listdir(bin_dir) if
                    os.path.isfile(os.path.join(bin_dir, _pyd)) and
                    os.path.splitext(_pyd)[0].startswith(PACKAGE_NAME) and
                    os.path.splitext(_pyd)[1] in [".pyd", ".so"]][0]
        shutil.move(pyd_path, os.path.join(os.path.split(pyd_path)[0], 'CFML_api', os.path.split(pyd_path)[1]))

        # Mark the libs for installation, adding them to
        # distribution.data_files seems to ensure that setuptools' record
        # writer appends them to installed-files.txt in the package's egg-info
        #
        # Also tried adding the libraries to the distribution.libraries list,
        # but that never seemed to add them to the installed-files.txt in the
        # egg-info, and the online recommendation seems to be adding libraries
        # into eager_resources in the call to setup(), which I think puts them
        # in data_files anyways.
        #
        # What is the best way?

        # These are the additional installation files that should be
        # included in the package, but are resultant of the cmake build
        # step; depending on the files that are generated from your cmake
        # build chain, you may need to modify the below code

        self.distribution.data_files = [os.path.join(self.install_dir,
                                                     os.path.basename(lib))
                                        for lib in libs]

        # Must be forced to run after adding the libs to data_files

        self.distribution.run_command("install_data")

        super().run()

        pass


class InstallCMakeScripts(install_scripts):
    """
    Install the scripts in the build dir
    """

    def run(self):
        """
        Copy the required directory to the build directory and super().run()
        """

        self.announce("Moving scripts files", level=3)

        # Scripts were already built in a previous step

        self.skip_build = True

        bin_dir = self.distribution.bin_dir
        scripts_dirs = []

        # scripts_dirs = [os.path.join(bin_dir, _dir) for _dir in
        #                 os.listdir(bin_dir) if
        #                 os.path.isdir(os.path.join(bin_dir, _dir))]
        #
        # for scripts_dir in scripts_dirs:
        #
        #     shutil.move(scripts_dir,
        #                 os.path.join(self.build_dir,
        #                              os.path.basename(scripts_dir)))

        # Mark the scripts for installation, adding them to
        # distribution.scripts seems to ensure that the setuptools' record
        # writer appends them to installed-files.txt in the package's egg-info

        self.distribution.scripts = scripts_dirs

        super().run()


class BuildCMakeExt(build_ext):
    """
    Builds using cmake instead of the python setuptools implicit build
    """

    def run(self):
        """
        Perform build_cmake before doing the 'normal' stuff
        """

        for extension in self.extensions:
            self.build_cmake(extension)
        # super().run()
        pass

    def build_cmake(self, extension: Extension):
        """
        The steps required to build the extension
        """
        cfg = "Debug" if self.debug else "Release"
        self.announce("Preparing the build environment", level=3)

        build_dir = pathlib.Path(self.build_temp)

        # don't store the abi3 info in filename
        extension._file_name = extension._full_name

        ext_path = pathlib.Path(self.get_ext_fullpath(extension.name))

        ext_ext = os.path.splitext(ext_path)[1]
        filename = extension._file_name + ext_ext
        extension_path = os.path.join(ext_path.parent, filename)

        os.makedirs(build_dir, exist_ok=True)
        os.makedirs(ext_path.parent.absolute(), exist_ok=True)

        # Now that the necessary directories are created, build

        self.announce("Configuring cmake project", level=3)

        # Change your cmake arguments below as necessary
        # Below is just an example set of arguments for building Blender as a Python module

        cmake_args = [
                         '-H' + SOURCE_DIR,
                         '-B' + self.build_temp,
                         f'-DPYTHON_INTERPRETER_PATH={sys.executable}',
                         python_lib_path,
                         f"-DCMAKE_Fortran_COMPILER={COMPILER}",
                         "-DARCH32=OFF",
                         "-DPYTHON_API=ON",
                         "-DUSE_HDF=OFF",
                         "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm
                         "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(os.path.join(build_dir), 'Release'),
                         "-DPYSETUP=ON"
                     ] + CFML_OVERRIDE
        build_args = []

        check_call(
            [get_cmake()] + cmake_args
        )

        self.announce("Building binaries", level=3)

        check_call(
            [get_cmake(), "--build", ".", '--target', 'install'] + build_args, cwd=self.build_temp
        )

        # Build finished, now copy the files into the copy directory
        # The copy directory is the parent directory of the extension (.pyd)
        self.announce("Moving built python module", level=3)

        bin_dir = os.path.join(os.getcwd(), COMPILER, 'Python_API', 'CFML_api')
        self.distribution.bin_dir = bin_dir

        pyd_path = [os.path.join(bin_dir, _pyd) for _pyd in
                    os.listdir(bin_dir) if
                    os.path.isfile(os.path.join(bin_dir, _pyd)) and
                    os.path.splitext(_pyd)[0].startswith(PACKAGE_NAME) and
                    os.path.splitext(_pyd)[1] in [".pyd", ".so"]][0]

        shutil.move(pyd_path, extension_path)
        # After build_ext is run, the following commands will run:
        #
        # install_lib
        # install_scripts
        #
        # These commands are subclassed above to avoid pitfalls that
        # setuptools tries to impose when installing these, as it usually
        # wants to build those libs and scripts as well or move them to a
        # different place. See comments above for additional information


class BuildCommand(distutils.command.build.build):

    def initialize_options(self):
        # this overrides the directory names for
        # build/lib and build/temp
        distutils.command.build.build.initialize_options(self)
        self.build_platlib = 'build/lib'
        self.build_temp = 'build/temp'


def build(kwargs: dict):
    ext_modules = [
        CMakeExtension(name=PACKAGE_NAME)
    ]
    cmdclass = {
        'build':           BuildCommand,
        'build_ext':       BuildCMakeExt,
        'install_data':    InstallCMakeLibsData,
        'install_lib':     InstallCMakeLibs,
        'install_scripts': InstallCMakeScripts,
        'bdist_wheel':     Bdist_wheel
    }
    kwargs.update({'ext_modules': ext_modules, 'cmdclass': cmdclass})


if __name__ == "__main__":
    data = read_configuration("./pyproject.toml")
    keys = ['name', 'version', 'description', 'classifiers']
    setup_kwargs = {
        key: data['project'][key] for key in keys
    }
    setup_kwargs['long_description'] = open(data['project']['readme'], 'r').read()
    setup_kwargs['long_description_content_type'] = 'text/markdown'
    setup_kwargs['author'] = ','.join([item['name'] for item in data['project']['authors']])
    setup_kwargs['author_email'] = data['project']['authors'][0]['email']
    setup_kwargs['license'] = data['project']['license']['text']
    # setup_kwargs['install_requires'] = data['project']['dependencies']
    build(setup_kwargs)
    setup(**setup_kwargs)
