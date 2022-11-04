from core.iFunction import IFunction
from solver.solverUtils import  getOffsetFunctionZ

class OffsetSolver:
    _surface: IFunction = None

    def __init__(self, function: IFunction):
        self._function = function

    def getOffsetFunctionZ(self, Pp, r):
        return getOffsetFunctionZ(self._function , Pp, r)
