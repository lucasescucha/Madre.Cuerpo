import math
import numpy as np

from  surfaces.iSurface import ISurface

TABS_HEIGHT_MARGIN = 0.1

X_VERSOR = np.array([1, 0, 0])
Y_VERSOR = np.array([0, 1, 0])
Z_VERSOR = np.array([0, 0, 1])

ZERO_VECTOR = ORIGIN = np.array([0, 0, 0])

SOLID_OPERATIONS_TOLERANCE = 0.1

def calculatePiecesData(referenceSurface, surfaceDimensions, 
            piecesConfiguration, axis):
        
    start = 0 if axis == "x" else (-surfaceDimensions.depth/2)
    end = surfaceDimensions.width/2 if axis == "x" \
        else (surfaceDimensions.depth/2) 

    lenght = referenceSurface.arcLenght(axis, start, end)

    pieceDimension = piecesConfiguration.width if axis == "x" \
        else piecesConfiguration.height

    pieces = math.ceil(lenght / pieceDimension)
    return start, end, pieces, lenght / pieces

def getShiftedClearPoygon(surface : ISurface, normal, dsplModule, 
        polygon : np.array, direction, clearence, fastVersion = True):
    
    displacement = normal * dsplModule * (1 if direction == "up" else -1)

    if fastVersion:
        return polygon + 1*displacement

    comparer = lambda polp, sfcp: (polp > (sfcp+clearence) if direction == "up" else  polp < (sfcp-clearence))

    while any([not comparer(v[2], surface.F(v[0:2])) for v in polygon]):
        polygon = polygon + displacement
    
    return polygon
