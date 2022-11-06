import numpy as np
import pyclipper as clipper

import surfaces.utils as sfcUtils
import utils.utils as utlUtils
import svg.utils as svgUtils

from solid.operations import getLeadMesh
from surfaces.iSurface import ISurface
from generator.utils import ORIGIN, Y_VERSOR, Z_VERSOR, ZERO_VECTOR, calculatePiecesData, getShiftedClearPoygon
from utils.utils import Configuration

def generateTabsShells(configuration: Configuration,
    referenceSurface : ISurface,  topOffsetSurface : ISurface, 
    baseOffsetSurface : ISurface):
    
    def getExtrudedPolygon(tabPolygon, bottomSurface, topSurface, thickness, n_u, clearence):
        bTabPolygon = getShiftedClearPoygon(bottomSurface, n_u, 
                thickness, tabPolygon, "down", clearence)
        
        tTabPolygon = getShiftedClearPoygon(topSurface, n_u, 
                thickness, tabPolygon, "up", clearence)
        
        return getLeadMesh(bTabPolygon, tTabPolygon, True) 

    def getTabBasePolygon(tabPolygon, flangeSize):
        FLOAT_INT_CONVERSION_FACTOR = 1000

        pco = clipper.PyclipperOffset()
        pco.AddPath(tabPolygon * FLOAT_INT_CONVERSION_FACTOR, 
                    clipper.JT_ROUND, clipper.ET_CLOSEDPOLYGON)
    
        flangeSize = flangeSize * FLOAT_INT_CONVERSION_FACTOR
        
        return np.array(pco.Execute(flangeSize)[0]) / FLOAT_INT_CONVERSION_FACTOR
        
    def locateTabPolygon(tabPolygon, lx, ly, xstart, ystart):
        referenceDirection = Y_VERSOR

        tabCenterPosition = np.array(
            [sfcUtils.arcLenght2Coordinate(referenceSurface, lx, "x", xstart),
             sfcUtils.arcLenght2Coordinate(referenceSurface, ly, "y", ystart)])

        n_u = sfcUtils.getUnitNormalVector(referenceSurface, tabCenterPosition)
        axis, angle = utlUtils.getRotationAngleAndAxis(np.array(Z_VERSOR), n_u)
        
        if not np.array_equal(axis, ZERO_VECTOR):
            tabPolygon = [utlUtils.rotateAroundAxis(np.append(v, [0]), axis, angle) for v in tabPolygon]
            referenceDirection = utlUtils.rotateAroundAxis(Y_VERSOR, axis, angle)

        direction = utlUtils.planeLineIntersection(ORIGIN, n_u, Y_VERSOR, Z_VERSOR)
        axis, angle = utlUtils.getRotationAngleAndAxis(referenceDirection, direction)

        if not np.array_equal(axis, ZERO_VECTOR):
            tabPolygon = [utlUtils.rotateAroundAxis(v, axis, angle) for v in tabPolygon]
        
        return tabPolygon + np.append(tabCenterPosition, [referenceSurface.F(tabCenterPosition)]), n_u

    surfDims = configuration.surface.dimensions
    moldConfig = configuration.manufacture.mold
    puzzleConfig = moldConfig.puzzle
    pDimensionsConfig = moldConfig.panels.dimensions

    pnlWidth, pnlHeight = pDimensionsConfig.width, pDimensionsConfig.height

    thickness = configuration.manufacture.mold.puzzle.thickness

    puzzleTabPolygon = svgUtils.getTabJoinPolygonFromSVGFile(
        puzzleConfig.tabs.file, puzzleConfig.tabs.width)

    puzzleTabBasePolygon = getTabBasePolygon(puzzleTabPolygon, 
        puzzleConfig.tabs.flangeSize)
    
    xstart, _, xPieces, pieceWidthLenght = calculatePiecesData(referenceSurface, 
        surfDims, puzzleConfig.pieces, "x")

    ystart, _, yPieces, pieceDepthLenght = calculatePiecesData(referenceSurface, 
        surfDims, puzzleConfig.pieces, "y")

    #Coloca los tab en el eje y
    for ix in range(xPieces):
        for iy in range(yPieces):

            if (iy+1) % pnlHeight == 0:
                continue

            lx = (ix+1/2)*pieceWidthLenght
            ly = (iy+1)*pieceDepthLenght

            tabPolygon, n_u = locateTabPolygon(puzzleTabPolygon, lx, ly, xstart, ystart)
            tabShell = getExtrudedPolygon(tabPolygon, referenceSurface, topOffsetSurface, thickness, n_u, thickness)

            tabPolygon, n_u = locateTabPolygon(puzzleTabBasePolygon, lx, ly, xstart, ystart)
            tabBaseShell = getExtrudedPolygon(tabPolygon, topOffsetSurface, baseOffsetSurface, thickness, n_u, thickness)

            yield tabShell, tabBaseShell

    R = np.array(((0, -1), (1,  0)))

    puzzleTabPolygon = np.dot(puzzleTabPolygon, R.T)
    puzzleTabBasePolygon = np.dot(puzzleTabBasePolygon, R.T)

    #Coloca los tab en el eje x
    for iy in range(yPieces):
        for ix in range(xPieces):

            if (ix+1) % pnlWidth == 0:
                continue

            ly = (iy+1/2)*pieceDepthLenght
            lx = (ix+1)*pieceWidthLenght

            tabPolygon, n_u = locateTabPolygon(puzzleTabPolygon, lx, ly, xstart, ystart)
            tabShell = getExtrudedPolygon(tabPolygon, referenceSurface, topOffsetSurface, thickness, n_u, thickness)

            tabPolygon, n_u = locateTabPolygon(puzzleTabBasePolygon, lx, ly, xstart, ystart)
            tabBaseShell = getExtrudedPolygon(tabPolygon, topOffsetSurface, baseOffsetSurface, thickness, n_u, thickness)

            yield tabShell, tabBaseShell
