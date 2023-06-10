import numpy as np

import surfaces.utils as sfcUtils
import FreeCADTools.utils as FreeCADUtils

from generator.utils import calculatePiecesData
from solid.operations import getLeadMesh
from surfaces.iSurface import ISurface
from surfaces.utils import calculateSamplingPoints
from utils.utils import Configuration

SAMPLING_POINTS_EXTRA_STEPS = 4
GRID_SIZE_FACTOR = 0.6
CUT_SURFACE_EXTRA_LENGHT = 0.05

def generatePiecesCutSurfaces(configuration: Configuration, 
    referenceSurface : ISurface):

    def calculateTopAndBottomPoints(x, y, thickness):
        position2d = [x, y]
        position3d = np.append(position2d, [referenceSurface.F(position2d)]) 
        
        n_u = sfcUtils.getUnitNormalVector(referenceSurface, position2d)

        return [position3d - (n_u*thickness*CUT_SURFACE_EXTRA_LENGHT), 
                position3d + (n_u*thickness*(1 + CUT_SURFACE_EXTRA_LENGHT))]

    surfDims = configuration.surface.dimensions
    puzzleConfig = configuration.manufacture.mold.puzzle
    grid = configuration.manufacture.grid

    thickness = configuration.manufacture.mold.puzzle.thickness

    xstart, xend, xPieces, pieceWidthLenght = calculatePiecesData(referenceSurface, 
        surfDims, puzzleConfig.pieces, "x")
    xVect = calculateSamplingPoints(xstart, xend, referenceSurface, "x", 
                grid.width*GRID_SIZE_FACTOR, SAMPLING_POINTS_EXTRA_STEPS)
    
    ystart, yend, yPieces, pieceDepthLenght = calculatePiecesData(referenceSurface, 
        surfDims, puzzleConfig.pieces, "y")
    yVect = calculateSamplingPoints(ystart, yend, referenceSurface, "y", 
                grid.height*GRID_SIZE_FACTOR, SAMPLING_POINTS_EXTRA_STEPS)

    for ix in range(1, xPieces):
        x = sfcUtils.arcLenght2Coordinate(referenceSurface, 
                    ix*pieceWidthLenght, "x", xstart)
        
        vertices = np.array(
                [calculateTopAndBottomPoints(x, y, thickness) for y in yVect])

        yield FreeCADUtils.convertMeshToSolid(getLeadMesh(vertices[:, 0], vertices[:, 1]))

    for iy in range(1, yPieces):
        y = sfcUtils.arcLenght2Coordinate(referenceSurface, iy*pieceDepthLenght, "y", ystart)
        
        vertices = np.array(
                [calculateTopAndBottomPoints(x, y, thickness) for x in xVect])

        yield FreeCADUtils.convertMeshToSolid(getLeadMesh(vertices[:, 0], vertices[:, 1]))

