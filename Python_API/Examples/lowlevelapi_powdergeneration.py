import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import numpy as np
import matplotlib.pyplot as plt

import CFML_api

filename = sys.argv[1]

#Read Cell, Space Group, Atom list and job info from CIF Files
cif_file = CFML_api.CIFFile(filename)

cell = cif_file.cell
space_group = cif_file.space_group
atom_list = cif_file.atom_list
job_info = cif_file.job_info

#Print the description of the different components
# cell.print_description()
# space_group.print_description()
# atom_list.print_description()
job_info.print_description()


#Possible to create/change an atom  form string (cf. Examples/test_atom_list.py)
# Create list from string
print("========\nCreate atom_list from string")
dat = [
'loop_                     ',
'_atom_site_label          ',
'_atom_site_fract_x        ',
'_atom_site_fract_y        ',
'_atom_site_fract_z        ',
'_atom_site_U_iso_or_equiv ',
'Sr 0.00000 0.00000 0.25000 0.00608',
'Ti 0.50000 0.00000 0.00000 0.00507',
'O1 0.00000 0.50000 0.25000 0.01646',
'O2 0.75000 0.25000 0.00000 0.02026']
atom_list = CFML_api.AtomList(dat)
atom_list.print_description()

# Switch atoms
atom1 = atom_list[0]
atom2 = atom_list[1]
atom_list[0]=atom2
atom_list[1]=atom1
print("========\nSwitched atom list")
atom_list.print_description()

# Modify atom position
atom = atom_list[0]
atom.xyz=np.asarray([0.2, 0.3, 0.4])
atom_list[0] = atom
print("========\nModify 1st atom position")
atom_list.print_description()

# Modify atom with string
atom = atom_list[0]
atom.from_string('ATOM C C 0.0 0.0 0.0 0.2 0.5')
atom_list[0] = atom
print("========\nModified 1st atom")
atom_list.print_description()

# Create atom from string
# ATOM   label    Chem/Scatt    x(s_x)   y(s_y)   z(s_z)   B_iso(s_b)    Occ(s_o)   mom  charge # String_with_info
atom2 = CFML_api.Atom('ATOM C C 0.2(0.1) 0.4(0.3) 0.5(0.4) 0.2 0.5(0.1)  3  4 #info')
print("========\nAtom created from string")
print(atom2)

# Replace atom
print("========\nReplace 4th atom in atom list")
atom_list[3] = atom2
atom_list.print_description()

#When creating or modifying the atom list from string (independently of the cell and space group
#there are a few subtleties concerning the multiplicities, occupancies and anisotropic dispacements

#Update the multiplicity and modify the occupancy when reading from a CIF array,
#in order to be in agreement with the definitions of Sfac in CrysFML
#This needs to be run when the AtomList is initialised from a CIF-like stringarray
#rather than from the class CIFFile
atom_list.set_mult_occ_cif(space_group)

#Update all Atom fields for Us, Bs and betas
#This needs to be run when the AtomList is initialised from a CIF-like stringarray
atom_list.set_all_anisotropic_displacement_parameters_cif(cell)

atom_list.print_description()
# Also possible to create job from string

print("========\nCreate job_info from string")
dat = [
'Title SrTiO3',
'Npatt 1',
'Patt_1 XRAY_2THE  1.54056    1.54056    1.00      0.0        135.0',
'UVWXY        0.025  -0.00020   0.01200   0.00150  0.00465',
'STEP         0.05 ',
'Backgd       50.000']


job_info_str = CFML_api.JobInfo(dat)
job_info_str.print_description()

job_info_str.range_2theta = (0.0, 12.0)

job_info_str.pattern_type = "NEUT_2THE"
job_info_str.lambda_ratio = 0.5

job_info_str.print_description()

print("========\n")
print("Compute Relection list\n")

# Friedel's pair F(h,k,l) = F(-h,-k,-l) in absence of anomalous dispersion phasing techniques
reflection_list = CFML_api.ReflectionList(cell, space_group, True, job_info)

reflection_list.compute_structure_factors(space_group, atom_list, job_info)

print(reflection_list.nref)

ref = reflection_list[0]

print(ref.hkl, ref.fcalc, ref.stl)

diffraction_pattern = CFML_api.DiffractionPattern(job_info,
                                                  reflection_list, cell.reciprocal_cell_vol)

#plt.plot(diffraction_pattern.x, diffraction_pattern.ycalc)
#plt.show()
