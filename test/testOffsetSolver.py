import unittest
import math
import numpy as np

from solver.solverUtils import getNormalHat, getNormalVector, solveF_Factory

class Test_TestOffsetSolver(unittest.TestCase):

    def test_getNormalVector(self):
        sqrt3 = math.sqrt(3)

        assert np.linalg.norm(
            getNormalVector(np.array([1, 1]))) == sqrt3
        assert np.linalg.norm(
            getNormalVector(np.array([-1, -1]))) == sqrt3
        assert np.linalg.norm(
            getNormalVector(np.array([1, -1]))) == sqrt3
        assert np.linalg.norm(
            getNormalVector(np.array([0, 0]))) == math.sqrt(1)

    def test_getNormalHat(self):
        assert np.linalg.norm(
           getNormalHat(np.array([2, -5]))) == 1

    def test_solveF(self):
        sqrt3 = math.sqrt(3)

        P = np.array([0, 0, 0])
        Pp = np.array([1/sqrt3, 1/sqrt3, 1/sqrt3])
        grad = np.array([-1, -1])

        solveF = solveF_Factory(Pp[0:2], lambda x: grad[0:2], 1)

        assert np.linalg.norm(solveF(P[0:2])) < 1e-5

    def test_svgpathtools(self):
        assert True

if __name__ == '__main__':
    unittest.main()
