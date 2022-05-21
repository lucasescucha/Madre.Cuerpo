import unittest
import FreeCADTools.utils as utils

class Test_TestFreeCADLoading(unittest.TestCase):

    def test_loadFreeCAD(self):
        doc = utils.createNewDocument()
        del doc

if __name__ == '__main__':
    unittest.main()
