import math

import numpy as np
import scipy.optimize as opt

from surfaces.iSurface import ISurface

def getUnitVector(vector):
    return vector / np.linalg.norm(vector)

def getUnitNormalVector(bottomSurface: ISurface, p):
    n = np.append(bottomSurface.gradF(p) * -1, 1)
    return getUnitVector(n)

def arcLenght2Coordinates(surface, lenghts, starts):
    x = arcLenght2Coordinate(surface, lenghts[0], "x", starts[0])
    y = arcLenght2Coordinate(surface, lenghts[1], "y", starts[1])

    return np.array([x, y])

def arcLenght2Coordinate(surface, lenght, axis, start):
    fun = lambda v: surface.arcLenght(axis, start, v[0]) - lenght
    return opt.root(fun, [0]).x[0]

def calculateSamplingPoints(astart: float, aend: float, surface : ISurface, 
            axis: str, gridAxisDimension: float, extraSteps = 0):
    
    curveLenght = surface.arcLenght(axis, astart, aend)
    pointsCount = math.ceil(curveLenght / gridAxisDimension)
    stepLenght = curveLenght / pointsCount 
    
    samplingPoints = []
    for i in range(-extraSteps, pointsCount + 1 + extraSteps):
        samplingPoints.append(arcLenght2Coordinate(surface, i * stepLenght, axis, astart))
    
    return np.array(samplingPoints)