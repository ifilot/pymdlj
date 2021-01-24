import unittest
from pymdlj import PyMDLJ
import numpy as np
import os

class TestSimulation(unittest.TestCase):

    def testSimulation(self):
        """
        Test simple simulation
        """
        pymdlj = PyMDLJ()
        results = pymdlj.simulate(os.path.join(os.path.dirname(__file__), 'settings', 'default.param.in'))

        self.assertTrue(len(results['ekin']) == results['nrsteps'])
        self.assertTrue(len(results['epot']) == results['nrsteps'])
        self.assertTrue(len(results['etot']) == results['nrsteps'])
        self.assertTrue(len(results['positions']) == results['nrparticles'])

if __name__ == '__main__':
    unittest.main()
