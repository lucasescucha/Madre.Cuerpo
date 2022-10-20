import unittest
import numpy as np

from utils.utils import getRotationAngleAndAxis
from utils.utils import rotateAroundAxis
from utils.utils import planeLineIntersection

class Test_TestUtils(unittest.TestCase):

    def test_getRotationAngleAndAxis(self):
        initialVector = [0, 0, 1] 
        directionVector = np.array([1, 0, 0])
         
        axis, angle = getRotationAngleAndAxis(initialVector, directionVector)
        
        finalVector = rotateAroundAxis(initialVector, axis, angle)

        assert np.isclose(directionVector, finalVector).all()

    def test_planeLineIntersection(self):
        intersection = planeLineIntersection([0, 0, 1], [0, 0, 1], [1, 0, 0], [0, 0, 1])
        assert np.isclose(intersection, [1, 0, 1]).all()

if __name__ == '__main__':
    unittest.main()
