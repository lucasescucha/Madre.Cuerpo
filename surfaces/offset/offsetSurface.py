from surfaces.iSurface import ISurface
from solver.offsetSurfaceSolver import OffsetSurfaceSolver
import numpy as np


class OffsetSurface(ISurface):
    def __init__(self, surface: ISurface, offset: float) -> None:
        self.solver = OffsetSurfaceSolver(surface)
        self.offset = offset

    def F(self, P: np.array) -> float:
        return self.solver.getOffsetSurfaceZ(P, self.offset)

    def gradF(self, P: np.array) -> np.array:
        raise NotImplementedError
