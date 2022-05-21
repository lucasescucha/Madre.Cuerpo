from surfaces.iSurface import ISurface

import scipy.optimize as opt
import numpy as np

SOLVER_TOLERANCE = 1e-5


class OffsetSurfaceSolver:
    _surface: ISurface = None

    def __init__(self, surface: ISurface):
        self._surface = surface

    @staticmethod
    def getNormalVector(gradF_P):
        return np.append((gradF_P * -1), [1])

    @staticmethod
    def getNormalHat(gradF_P):
        normal = OffsetSurfaceSolver.getNormalVector(gradF_P)
        return normal/np.linalg.norm(normal)

    @staticmethod
    def solveF(P, Pp, gradF, r):
        n_hat = OffsetSurfaceSolver.getNormalHat(gradF(P))
        return Pp - P - (n_hat[0:2] * r)

    def getOffsetSurfaceZ(self, Pp, r):
        surfaceGradF = self._surface.gradF

        def solveFWrapper(P): return OffsetSurfaceSolver.solveF(
            P, Pp, surfaceGradF, r)
        Pr = opt.excitingmixing(solveFWrapper, Pp, f_tol=SOLVER_TOLERANCE)

        n_hat = OffsetSurfaceSolver.getNormalHat(surfaceGradF(Pr))
        PpZ = self._surface.F(Pr) + (n_hat[2] * r)

        return PpZ
