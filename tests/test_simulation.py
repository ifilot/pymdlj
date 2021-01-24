import unittest
from pymdlj import PyMDLJ
import numpy as np
import os

class TestSimulation(unittest.TestCase):

    def testSimulation(self):
        """
        Test Isosurface Generation of a Gaussian
        """
        pymdlj = PyMDLJ()
        pymdlj.simulate(os.path.join(os.path.dirname(__file__), 'settings', 'default.param.in'))

if __name__ == '__main__':
    unittest.main()
