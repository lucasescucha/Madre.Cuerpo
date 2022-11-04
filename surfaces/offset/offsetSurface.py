from surfaces.iSurface import ISurface
from solver.offsetSolver import OffsetSolver
import numpy as np


class OffsetSurface(ISurface):
    def __init__(self, surface: ISurface, offset: float) -> None:
        self.solver = OffsetSolver(surface)
        self.offset = offset

    def F(self, P: np.array) -> float:
        return self.solver.getOffsetFunctionZ(P, self.offset)

    def gradF(self, P: np.array) -> np.array:
        raise NotImplementedError

    def arcLenght(self, direction: str, start: float, end: float) -> float:
        raise NotImplementedError
