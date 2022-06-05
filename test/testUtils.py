import unittest
import numpy as np

from utils.utils import getRotationAngleAndAxis
from utils.utils import rotateAroundAxis

class Test_TestUtils(unittest.TestCase):

    def test_getRotationAngleAndAxis(self):
        initialVector = [0, 0, 1] 
        directionVector = np.array([1, 0, 0])
         
        axis, angle = getRotationAngleAndAxis(initialVector, directionVector)
        
        finalVector = rotateAroundAxis(initialVector, axis, angle)

        assert np.isclose(directionVector, finalVector).all()

if __name__ == '__main__':
    unittest.main()
