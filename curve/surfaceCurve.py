from curve.iCurve import ICurve
from solver.solverUtils import getNormalHat
from surfaces.iSurface import ISurface

import numpy as np

class SurfaceCurve(ICurve):
    def __init__(self, surface: ISurface, axis:str, secAxisValue:float) -> None:
        self._surface = surface
        self._axis = axis
        self._secAxisValue = secAxisValue

    @staticmethod
    def getSurfacePoint(axis, P, secAxisValue):
        return [P[0], secAxisValue] if axis=="x" else [secAxisValue, P[0]]

    def F(self, P: np.array) -> float:
        Pf = SurfaceCurve.getSurfacePoint(self._axis, P, self._secAxisValue)
        return self._surface.F(Pf)

    def gradF(self, P: np.array) -> np.array:
        mainAxisIndex = (0 if self._axis=="x" else 1)
        Pf = SurfaceCurve.getSurfacePoint(self._axis, P, self._secAxisValue)
        return np.array([self._surface.gradF(Pf)[mainAxisIndex]])
    
    def arcLenght(self, start: float, end: float) -> float:
        raise NotImplementedError

    def FOffset(self, P: np.array, r: float) -> np.array:
        n_hat = getNormalHat(self.gradF(P))
        return np.append(P, [self.F(P)]) + n_hat*r