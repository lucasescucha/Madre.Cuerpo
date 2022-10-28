import unittest

from surfaces.meshSurface import checkPointInTriangle

class Test_TestMeshSurfacesUtils(unittest.TestCase):

    def test_checkPointInTriangle(self):
        triangle = [[-1, 0], [1, 1], [-1, -1]]

        P = [0, 0]
        assert checkPointInTriangle(P, triangle) == True

        P = [2, 2]
        assert checkPointInTriangle(P, triangle) == False

if __name__ == '__main__':
    unittest.main()
