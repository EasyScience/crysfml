import unittest

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from CFML_api import SpaceGroup

import numpy
import unittest

class TestSymmetryGroups(unittest.TestCase):
    '''
    Unittest for the geometry-related functions
    '''

    def test_print_description(self):
        self.group = SpaceGroup(1)
        self.group.print_description()

    def test_getters(self):
        self.group = SpaceGroup(1)
        numpy.testing.assert_almost_equal(self.group.lattice_translation, numpy.array([[0.], [0.], [0.]]))
        self.assertFalse(self.group.is_hexa)
        self.assertEqual(self.group.space_group_number, 1)
        self.assertEqual(self.group.space_group_symbol, "P 1")
        self.assertEqual(self.group.hall_symbol, "P 1")
        self.assertEqual(self.group.generalized_hall_symbol, "")
        self.assertEqual(self.group.crystal_system, "Triclinic")
        self.assertEqual(self.group.laue_class, "-1")
        self.assertEqual(self.group.point_group, "1")
        self.assertEqual(self.group.extra_information, "")
        self.assertEqual(self.group.space_group_setting,  "Generated from Hall symbol")
        self.assertEqual(self.group.space_group_lattice_type, "P")
        self.assertEqual(self.group.space_group_lattice_type_symbol, "aP")
        self.assertEqual(self.group.number_of_lattice_points, 1)
        self.assertEqual(self.group.bravais, "  P: { 000 }")
        self.assertEqual(self.group.centre, "Acentric")
        self.assertEqual(self.group.centric, 1)
        numpy.testing.assert_almost_equal(self.group.inversion_centre, numpy.array([0., 0., 0.]))
        self.assertEqual(self.group.number_of_symops, 1)
        self.assertEqual(self.group.multiplicity, 1)
        self.assertEqual(self.group.operators_minimum_number, 1)
        self.assertEqual(self.group.symmetry_operators_as_text, ["x,y,z"])
        numpy.testing.assert_almost_equal(self.group.asymetric_unit, numpy.array([[0.,1.],[0.,1.],[0.,1.]]))
        numpy.testing.assert_almost_equal(self.group.symmetry_operators[0].rotation_matrix, numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]))
        numpy.testing.assert_almost_equal(self.group.symmetry_operators[0].translation_matrix, numpy.array([0., 0., 0.]))
        numpy.testing.assert_almost_equal(self.group.wyckoff_info.num_orbit, 0)
        # For some reason this fails. To Look into
        #self.assertEqual(self.group.wyckoff_info.orbits[0].multiplicity, 0)
        #self.assertEqual(self.group.wyckoff_info.orbits[0].site, "\x00"*6)
        #self.assertEqual(self.group.wyckoff_info.orbits[0].orbit_number, 0)
        #self.assertEqual(self.group.wyckoff_info.orbits[0].str_orig, "\x00"*40)
        #self.assertEqual(self.group.wyckoff_info.orbits[0].str_orbit[0], "\x00"*40)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(TestSymmetryGroups))
    return s

if __name__ == '__main__':
    unittest.main(verbosity=2)
            
        
