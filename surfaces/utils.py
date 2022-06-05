import numpy as np

from surfaces.iSurface import ISurface

def getUnitNormalVector(bottomSurface: ISurface, p):
    n = np.append(bottomSurface.gradF(p) * -1, 1)
    return n / np.linalg.norm(n)