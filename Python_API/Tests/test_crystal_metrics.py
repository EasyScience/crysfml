import unittest

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from CFML_api import Cell

import numpy as np
import unittest

class TestCell(unittest.TestCase):
    ''' 
    Unittest for the crystal cell class
    '''

    def test_print_description(self):
        abc  = np.asarray([5.,5.,8.],dtype='float32')
        angl = np.asarray([90.,90.,80.],dtype='float32')
        self.cell = Cell(abc,angl)
        self.cell.print_description()

    def test_getters(self):
        abc  = np.asarray([5.,5.,8.],dtype='float32')
        angl = np.asarray([90.,90.,80.],dtype='float32')
        self.cell = Cell(abc,angl)

        # Relax the comparison to 4dp, MuscleLinux complains :-( 

        np.testing.assert_almost_equal(self.cell.lattpar, np.array([5., 5., 8.]))
        np.testing.assert_almost_equal(self.cell.lattangle, np.array([90., 90., 80.]), decimal=4)
        np.testing.assert_almost_equal(self.cell.lsq_lattpar, np.array([0, 0, 0]), decimal=4)
        np.testing.assert_almost_equal(self.cell.lsq_lattangle, np.array([0, 0, 0]), decimal=4)
        np.testing.assert_almost_equal(self.cell.lattpar_std_dev, np.array([0, 0, 0]), decimal=4)
        np.testing.assert_almost_equal(self.cell.lattangle_std_dev, np.array([0, 0, 0]), decimal=4)
        np.testing.assert_almost_equal(self.cell.reciprocal_lattpar, np.array([0.2030853, 0.2030853, 0.125000 ]), decimal=4)
        np.testing.assert_almost_equal(self.cell.reciprocal_lattangle, np.array([90., 90., 100.]), decimal=4)
        np.testing.assert_almost_equal(self.cell.direct_metric_tensor, np.array([[25., 4.3412046, 0.],[4.3412046, 25., 0.],[0.,0.,64.]]), decimal=4)
        np.testing.assert_almost_equal(self.cell.reciprocal_metric_tensor,
                                       np.array([[0.04124364, -0.0071618841 , 0.],[-0.0071618841 , 0.04124364, 0.],[0.,0.,0.015625]]), 
                                       decimal=4)
        np.testing.assert_almost_equal(self.cell.crystal_to_orth_matrix,
                                       np.array([[4.924039, 0., 0.],[0.8682409, 5., 0.],[0.,0.,8.]]), decimal=4)
        np.testing.assert_almost_equal(self.cell.orth_to_crystal_matrix,
                                       np.array([[0.20308532, 0., 0.],[-0.035265397,  0.2, 0.],[0., 0., 0.125]]), decimal=4)
        
                                       
def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(TestCell))
    return s

if __name__ == '__main__':
    unittest.main(verbosity=2)
