import math
import operator
import os
import shutil

import numpy as np
import pyclipper as clipper
from solid.operations import getLeadMesh, getShellMesh, surfaceToMeshSurface, surfaceToOffsetMeshSurface

from surfaces.meshSurface import MeshSurface

import svg.utils as svgUtils
import surfaces.utils as sfcUtils
import utils.utils as utlUtils
import FreeCADTools.utils as FreeCADUtils

from configuration import configurationData
from surfaces.iSurface import ISurface

from surfaces.motherSurface import MotherSurface
from surfaces.shiftDecorator import ShiftDecorator
from surfaces.utils import calculateSamplingPoints

from utils.utils import Configuration, rotateAroundAxis

TABS_HEIGHT_MARGIN = 0.1

X_VERSOR = [1, 0, 0]
Y_VERSOR = [0, 1, 0]
Z_VERSOR = [0, 0, 1]

ZERO_VECTOR = ORIGIN = [0, 0, 0]

SOLID_OPERATIONS_TOLERANCE = 0.1

def generateMeshSurfaces(configuration: Configuration):
    def calculateReferenceSurface(surfaceConfig):
        zOffset = surfaceConfig.shape.zeroHeight

        # z = zOffset + a*x^2 -  b*y^2

        # y=0 -> z = zOffset + a*x^2 -> height = zOffset + a*(width/2)^2
        # a = (height - zOffset) / (width/2)^2
        a = (surfaceConfig.dimensions.height - zOffset) / \
            ((surfaceConfig.dimensions.width/2)**2)

        # x=0 -> z = zOffset - b*y^2 -> 0 = zOffset - b*(depth/2)^2
        # b = zOffset / (depth/2)^2
        b = zOffset / ((surfaceConfig.dimensions.depth/2)**2)

        return ShiftDecorator(MotherSurface(a, b),  np.array([0, 0, zOffset]))

    surface = configuration.surface.dimensions
    grid = configuration.manufacture.grid

    thickness = configuration.manufacture.mold.puzzle.thickness
    tabsThickness = configuration.manufacture.mold.puzzle.tabs.thickness

    referenceSurface = calculateReferenceSurface(configuration)

    xstart, xend = 0, surface.width/2
    xVect = calculateSamplingPoints(xstart, xend, referenceSurface, "x", grid.width)

    ystart, yend = (-surface.depth/2), (surface.depth / 2)
    yVect = calculateSamplingPoints(ystart, yend, referenceSurface, "y", grid.height)

    bottomSurface = surfaceToMeshSurface(referenceSurface, xVect, yVect)

    topSurface = surfaceToOffsetMeshSurface(bottomSurface, thickness, xVect, yVect)

    baseOffset = thickness + tabsThickness
    baseSurface = surfaceToOffsetMeshSurface(bottomSurface, baseOffset, xVect, yVect)

    return referenceSurface, bottomSurface, topSurface, baseSurface

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

def getShiftedClearPoygon(surface : MeshSurface, normal, dsplModule, 
        polygon : np.array, direction, clearence, fastVersion = True):
    
    displacement = normal * dsplModule * (1 if direction == "up" else -1)

    if fastVersion:
        return polygon + 10*displacement

    comparer = lambda polp, sfcp: (polp > (sfcp+clearence) if direction == "up" else  polp < (sfcp-clearence))

    while any([not comparer(v[2], surface.F(v[0:2])) for v in polygon]):
        polygon = polygon + displacement
    
    return polygon

def generateTabsShells(configuration: Configuration,
    referenceSurface : ISurface,  bottomSurface: MeshSurface, 
    topSurface : MeshSurface, baseSurface : MeshSurface):
    
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
            tabShell = getExtrudedPolygon(tabPolygon, bottomSurface, topSurface, 3*thickness, n_u, thickness)

            tabPolygon, n_u = locateTabPolygon(puzzleTabBasePolygon, lx, ly, xstart, ystart)
            tabBaseShell = getExtrudedPolygon(tabPolygon, topSurface, baseSurface, 3*thickness, n_u, thickness)

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
            tabShell = getExtrudedPolygon(tabPolygon, bottomSurface, topSurface, 3*thickness, n_u, thickness)

            tabPolygon, n_u = locateTabPolygon(puzzleTabBasePolygon, lx, ly, xstart, ystart)
            tabBaseShell = getExtrudedPolygon(tabPolygon, topSurface, baseSurface, 3*thickness, n_u, thickness)

            yield tabShell, tabBaseShell

def generatePiecesCutSurfaces(configuration: Configuration, 
    referenceSurface : ISurface,  bottomSurface: MeshSurface, 
    topSurface : MeshSurface):

    def calculateTopAndBottomPoints(x, y, thickness):
        position2d = [x, y]
        position3d = np.append(position2d, [referenceSurface.F(position2d)]) 
        
        n_u = sfcUtils.getUnitNormalVector(referenceSurface, position2d)
        
        bottomOffset = getShiftedClearPoygon(bottomSurface, n_u, 3*thickness, 
                [position3d], "down", thickness)

        topOffset = getShiftedClearPoygon(topSurface, n_u, 3*thickness, 
                [position3d], "up", thickness)
            
        return [bottomOffset, topOffset] 

    surfDims = configuration.surface.dimensions
    puzzleConfig = configuration.manufacture.mold.puzzle
    grid = configuration.manufacture.grid

    thickness = configuration.manufacture.mold.puzzle.thickness

    xstart, xend, xPieces, pieceWidthLenght = calculatePiecesData(referenceSurface, 
        surfDims, puzzleConfig.pieces, "x")
    xVect = calculateSamplingPoints(xstart, xend, referenceSurface, "x", grid.width)
    
    ystart, yend, yPieces, pieceDepthLenght = calculatePiecesData(referenceSurface, 
        surfDims, puzzleConfig.pieces, "y")
    yVect = calculateSamplingPoints(ystart, yend, referenceSurface, "y", grid.height)

    for ix in range(1, xPieces):
        x = sfcUtils.arcLenght2Coordinate(referenceSurface, ix*pieceWidthLenght, "x", xstart)
        
        vertices = np.array(
            [calculateTopAndBottomPoints(x, y, thickness) for y in yVect])

        yield getLeadMesh(vertices[:, 0], vertices[:, 1])

    for iy in range(1, yPieces):
        y = sfcUtils.arcLenght2Coordinate(referenceSurface, iy*pieceDepthLenght, "y", ystart)
        
        vertices = np.array(
            [calculateTopAndBottomPoints(x, y, thickness) for x in xVect])

        yield getLeadMesh(vertices[:, 0], vertices[:, 1])

def generatePanelsLeads(configuration, referenceSurface):
    LEFT_LEAD, RIGHT_LEAD, TOP_LEAD, BOTTOM_LEAD = 0, 1, 2, 3

    TOP_FLANGE, BOTTOM_FLANGE, LEFT_FLANGE, RIGHT_FLANGE = 0, 1, 2, 3

    MIDDLE_FLANGE, CORNER_FLANGE = 0, 1

    TOP_DIRECTION, RIGHT_DIRECTION = 1
    BOTTOM_DIRECTION, LEFT_DIRECTION = -1

    surfDims = configuration.surface.dimensions
    moldConfig = configuration.manufacture.mold
    puzzleConfig = moldConfig.puzzle

    moldThickness = moldConfig.puzzle.thickness

    leadHeight = moldConfig.panels.leads.height
    leadThickness = moldConfig.panels.leads.thickness
    flangeSize = moldConfig.panels.leads.flangeSize
    leadScrews = moldConfig.panels.leads.screws
    leadScrewsDiameter = moldConfig.panels.leads.screwsDiameter
    pieceInsertNutDiameter = moldConfig.panels.leads.insertNutDiameter
    pieceInsertNutDepth = moldConfig.panels.leads.insertNutDepth

    grid = configuration.manufacture.grid

    xstart, _, piecesX, pieceWidthLenght = calculatePiecesData(
        referenceSurface, surfDims, puzzleConfig.pieces, "x")
    
    ystart, _, piecesY, pieceDepthLenght = calculatePiecesData(
        referenceSurface, surfDims, puzzleConfig.pieces, "y")

    def getFlangeDirections(leadSide, flangePosition, flangeType, p0):
        n_u = sfcUtils.getUnitNormalVector(referenceSurface, p0)
        
        if leadSide == LEFT_LEAD:
            if flangePosition == TOP_FLANGE:
                t = getTangentVector(p0, "y", TOP_DIRECTION)
            elif flangePosition == BOTTOM_FLANGE:
                t = getTangentVector(p0, "x", LEFT_DIRECTION)
        
        if leadSide == RIGHT_LEAD:
            if flangePosition == TOP_FLANGE:
                t = getTangentVector(p0, "x", RIGHT_DIRECTION)
            elif flangePosition == BOTTOM_FLANGE:
                t = getTangentVector(p0, "y", BOTTOM_DIRECTION)
              
        if leadSide == TOP_LEAD:
            if flangePosition == LEFT_FLANGE:
                t = getTangentVector(p0, "y", TOP_DIRECTION)
            elif flangePosition == RIGHT_FLANGE:
                t = getTangentVector(p0, "x", RIGHT_DIRECTION)
        
        if leadSide == BOTTOM_LEAD:
            if flangePosition == LEFT_FLANGE:
                t = getTangentVector(p0, "x", LEFT_DIRECTION)
            elif flangePosition == RIGHT_FLANGE:
                t = getTangentVector(p0, "y", BOTTOM_DIRECTION)

        offsetNedeed = \
            (leadSide == LEFT_LEAD and flangePosition == TOP_FLANGE) or \
            (leadSide == RIGHT_LEAD and flangePosition == BOTTOM_FLANGE) or \
            (leadSide == TOP_LEAD and flangePosition == RIGHT_FLANGE) or \
            (leadSide == BOTTOM_LEAD and flangePosition == LEFT_FLANGE)
        
        if flangeType == MIDDLE_FLANGE:
            d1 = rotateAroundAxis(t, n_u, np.pi/2) if offsetNedeed else t
        elif flangeType == CORNER_FLANGE:
            d1 = rotateAroundAxis(t, n_u, np.pi/4)

        d2 = rotateAroundAxis(d1, n_u, np.pi/2 if offsetNedeed else -np.pi/2)

        return d1, d2

    def getOffsetPoints(points, r):
        return [referenceSurface.FOffset(p, r) for p in points]

    def getTangentVector(p, axis, direction):
        grad = np.abs(referenceSurface.gradF(p)) * direction
        t = np.array([1, 0, grad[0]] if axis == "x" else [0, 1, grad[1]])
        t_hat = t/np.linalg.norm(t)
        
        return np.append(t_hat, [0])

    def extendVertices(vertices, axis, direction, distance):
        return [(v + (getTangentVector(v, axis, direction)*distance)) 
                for v in vertices]

    def getLeadPieceScrewsDrills(mainAxis, malstart, malend, salStart, 
                        saDir, isPiece = True):

        secondaryAxis = "y" if mainAxis == "x" else "x"

        pa_start = xstart if mainAxis == "x" else ystart
        sa_start = xstart if secondaryAxis == "x" else ystart

        sa_pos = sfcUtils.arcLenght2Coordinate(
                referenceSurface, salStart, secondaryAxis, sa_start)

        ldrills = (malend-malstart)/(leadScrews+1)
 
        for i in range(leadScrews):
            ldrill = malstart + (i+1)*ldrills 
            
            ma_pos = sfcUtils.arcLenght2Coordinate(
                referenceSurface, ldrill, mainAxis, pa_start)

            pos = [ma_pos, sa_pos] if mainAxis == "x" else [sa_pos, ma_pos]

            drillDirection = getTangentVector(pos, secondaryAxis, saDir)
            drillCenter = getOffsetPoints([pos], moldThickness - (leadHeight/2))

            if isPiece:
                drillCenter = drillCenter - drillDirection * pieceInsertNutDepth
                drillDepth = 2 * pieceInsertNutDepth
                drillRadius = pieceInsertNutDiameter/2
            else:
                drillCenter = drillCenter - drillDirection * leadThickness
                drillDepth = 3 * leadThickness
                drillRadius = leadScrewsDiameter/2

            yield FreeCADUtils.createCylinder(drillRadius, drillDepth, drillCenter, drillDirection)

    def getFlangeScrewDirll(leadSide, flangePosition, flangeType, lx, ly):
        drillDepth = 3*leadThickness
        drillRadius = leadScrewsDiameter/2

        p0 = sfcUtils.arcLenght2Coordinates(referenceSurface, [lx, ly], [xstart, ystart])

        v0 = getOffsetPoints([p0], moldThickness - (leadHeight/2))

        d1, d2 = getFlangeDirections(leadSide, flangePosition, flangeType, p0)

        drillCenter = v0 + (d1 * flangeSize/2) - (d2 * leadThickness)
        drillDirection = d2 

        return FreeCADUtils.createCylinder(drillRadius, drillDepth, drillCenter, drillDirection)

    def getLeadMeshAndDrills(mainAxis, malStart, malEnd, salStart, saDir):
        secondaryAxis = "y" if mainAxis == "x" else "x"

        pa_start = xstart if mainAxis == "x" else ystart
        sa_start = xstart if secondaryAxis == "x" else ystart
        
        gridDim = grid.width if mainAxis == "x" else grid.height

        mpa_start = sfcUtils.arcLenght2Coordinate(
                referenceSurface, malStart, mainAxis, pa_start)
        mpa_end = sfcUtils.arcLenght2Coordinate(
                referenceSurface, malEnd, mainAxis, pa_start)
        
        secondaryAxisPoint = sfcUtils.arcLenght2Coordinate(
                referenceSurface, salStart, secondaryAxis, sa_start)

        mainAxisVect = calculateSamplingPoints(mpa_start, mpa_end, 
                            referenceSurface, mainAxis, gridDim)

        if mainAxis == "x":
            points = [[x, secondaryAxisPoint] for x in mainAxisVect]
        else:
            points = [[secondaryAxisPoint, y] for y in mainAxisVect]

        rt = moldThickness
        rb = moldThickness - leadHeight

        topVertices = getOffsetPoints(rt, points)
        bottomVertices = getOffsetPoints(rb, points)

        leadAFace = [topVertices, bottomVertices]
        
        eTopVertices = extendVertices(topVertices, secondaryAxis, saDir, leadThickness)
        eBottomVertices = extendVertices(bottomVertices, secondaryAxis, saDir, leadThickness)

        leadBFace = [eTopVertices, eBottomVertices]

        drills = getLeadPieceScrewsDrills(mainAxis, malStart, malEnd, salStart, saDir, isPiece=False)

        return [utlUtils.getShellMesh(leadAFace, leadBFace), drills]

    def getFlangeMeshAndDrill(leadSide, flangePosition, flangeType, lx, ly):
        p0 = sfcUtils.arcLenght2Coordinates(referenceSurface, [lx, ly], [xstart, ystart])

        rt = moldThickness
        rb = moldThickness - leadHeight

        d1, d2 = getFlangeDirections(leadSide, flangePosition, flangeType, p0)

        v0a = [getOffsetPoints([p0], rb), getOffsetPoints([p0], rt)]        
        v0b = v0a + d1 * flangeSize
        
        v1a = v0a + d2 * leadThickness
        v1b = v0b + d2 * leadThickness
        
        drill = getFlangeScrewDirll(leadSide, flangePosition, flangeType, lx, ly)

        return [utlUtils.getShellMesh([v0a, v0b], [v1a, v1b]), [drill]]

    pnWidth = moldConfig.panels.dimensions.width
    pnHeight = moldConfig.panels.dimensions.height

    for ix in range(piecesX):        
        for iy in range(piecesY):
            rightLead, leftLead = [], []
            topLead, bottomLead = [], []

            piecesDrills = []

            leftSide = (ix % pnWidth) == 0
            rightSide = ((ix+1) % pnWidth) == 0

            topSide = (iy % pnHeight) == 0
            bottomSide = ((iy+1) % pnHeight) == 0

            # Side leads
            
            lystart = iy*pieceDepthLenght
            lyend = lystart + pieceDepthLenght

            if leftSide:               
                lxstart = ix*pieceWidthLenght
                
                leftLead.append(
                    getLeadMeshAndDrills("y", lystart, lyend, lxstart, LEFT_DIRECTION))  

                piecesDrills.append(
                    getLeadPieceScrewsDrills("y", lystart, lyend, lxstart, LEFT_DIRECTION))     
            else:
                if topSide:
                    #agregar oreja a 90° arriba a la izquierda sólo a top 
                    topLead.append(
                        getFlangeMeshAndDrill(TOP_LEAD, LEFT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
                
                if bottomLead:
                    #agregar oreja a 90° abajo a la izquierda sólo a bottom 
                    bottomLead.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, LEFT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
            
            if rightSide:
                lxstart = (ix+1)*pieceWidthLenght

                rightLead.append(
                    getLeadMeshAndDrills("y", lystart, lyend, lxstart, RIGHT_DIRECTION))
                
                piecesDrills.append(
                    getLeadPieceScrewsDrills("y", lystart, lyend, lxstart, RIGHT_DIRECTION))     
            else:
                if topSide:
                    #agregar oreja a 90° arriba a la derecha sólo a top 
                    topLead.append(
                        getFlangeMeshAndDrill(TOP_LEAD, RIGHT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
                
                if bottomLead:
                    #agregar oreja a 90° abajo a la derecha sólo a bottom 
                    bottomLead.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, RIGHT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))

            lxstart = ix*pieceWidthLenght
            lxend = lxstart + pieceWidthLenght

            if topSide:
                lystart = iy*pieceDepthLenght

                topLead.append(
                    getLeadMeshAndDrills("x", lxstart, lxend, lystart, TOP_DIRECTION))
                
                piecesDrills.append(
                    getLeadPieceScrewsDrills("x", lxstart, lxend, lystart, TOP_DIRECTION))     
            else:
                if leftSide:
                    #agregar oreja a 90° abajo a la izquierda sólo a left
                    leftLead.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, LEFT_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 
                
                if rightSide:
                    #agregar oreja a 90° abajo a la derecha sólo a right 
                    rightLead.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, RIGHT_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 

            if bottomSide:
                lystart = (iy+1)*pieceDepthLenght

                bottomLead.append(
                    getLeadMeshAndDrills("x", lxstart, lxend, lystart, BOTTOM_DIRECTION))

                piecesDrills.append(
                    getLeadPieceScrewsDrills("x", lxstart, lxend, lystart, BOTTOM_DIRECTION))
            else:
                if leftSide:
                    #agregar oreja a 90° arriba a la izquierda sólo a left 
                    leftLead.append(
                        getFlangeMeshAndDrill(TOP_LEAD, LEFT_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 
                
                if rightSide:
                    #agregar oreja a 90° arriba a la derecha sólo a right 
                    rightLead.append(
                        getFlangeMeshAndDrill(TOP_LEAD, RIGHT_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 

            if leftSide and topSide:
                #agregar oreja a 45° arriba a la izquierda a left y a top
                leftLead.append(
                        getFlangeMeshAndDrill(LEFT_LEAD, TOP_FLANGE, CORNER_FLANGE, lxstart, lystart))
                topLead.append(
                        getFlangeMeshAndDrill(TOP_LEAD, LEFT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 
            
            if leftSide and bottomSide:
                #agregar oreja a 45° abajo a la izquierda a left y a bottom
                leftLead.append(
                        getFlangeMeshAndDrill(LEFT_LEAD, BOTTOM_FLANGE, CORNER_FLANGE, lxstart, lystart))
                bottomLead.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, LEFT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            if rightSide and topSide:
                #agregar oreja a 45° arriba a la derecha a right y a top
                rightLead.append(
                        getFlangeMeshAndDrill(RIGHT_LEAD, TOP_FLANGE, CORNER_FLANGE, lxstart, lystart))
                topLead.append(
                        getFlangeMeshAndDrill(TOP_LEAD, RIGHT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            if rightSide and bottomSide:
                #agregar oreja a 45° abajo a la derecha a right y a bottom
                rightLead.append(
                        getFlangeMeshAndDrill(RIGHT_LEAD, BOTTOM_FLANGE, CORNER_FLANGE, lxstart, lystart))
                bottomLead.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, RIGHT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            yield rightLead, leftLead, topLead, bottomLead, piecesDrills

def checkParts(parts):
    for part in parts:
        if not part.isValid():
            return False

    return True

def run():
    configuration = Configuration(configurationData)

    if not configuration.check():
        raise ArithmeticError

    referenceSurface, bottomSurface, topSurface, baseSurface = \
        generateMeshSurfaces(configuration)

    surfaceTriangleMesh = getShellMesh(bottomSurface, topSurface)
    baseTriangleMesh = getShellMesh(topSurface, baseSurface)

    document = FreeCADUtils.createNewDocument()

    surfaceSolid = FreeCADUtils.convertMeshToSolid(surfaceTriangleMesh)
    baseSolid = FreeCADUtils.convertMeshToSolid(baseTriangleMesh)

    tabsShells = list(generateTabsShells(configuration, referenceSurface,  bottomSurface, 
                topSurface, baseSurface))
    
    tabsSolid = []

    body = surfaceSolid

    for shell in tabsShells:
        result = FreeCADUtils.slicePart(body, FreeCADUtils.createMesh(shell[0]))
        
        volumeSortedResult = sorted(
            result, key=operator.attrgetter("Volume"), reverse=True)
        
        body, tab = volumeSortedResult[0], volumeSortedResult[1:]

        result = FreeCADUtils.slicePart(baseSolid, FreeCADUtils.createMesh(shell[1]))

        volumeSortedResult = sorted(
            result, key=operator.attrgetter("Volume"), reverse=True)

        _, baseTab = volumeSortedResult[0], volumeSortedResult[1:]

        tabsSolid.append(tab.fuse(baseTab))

    def getSolid(leadData):
        leadSolid = None
        for part in leadData:
            solid = FreeCADUtils.convertMeshToSolid(part[0])
            solid = solid.cut(part[1])

            if leadSolid==None:
                leadSolid = solid
            else:
                leadSolid = leadSolid.fuse(solid)
        
        return leadSolid

    leadsSolids = []

    for leadsAndDrills in generatePanelsLeads(configuration, referenceSurface):
        rightLead, leftLead, topLead, bottomLead, surfaceDrills = leadsAndDrills

        if rightLead != []: 
            leadsSolids.append(getSolid(rightLead))
        if leftLead != []: 
            leadsSolids.append(getSolid(leftLead))
        if topLead != []: 
            leadsSolids.append(getSolid(topLead))
        if bottomLead != []: 
            leadsSolids.append(getSolid(bottomLead))
        
        for drill in surfaceDrills:
            body = body.cut(drill)

    cutSurfaces = generatePiecesCutSurfaces(configuration, referenceSurface, bottomSurface, topSurface)

    cutSurfacesMesh = [FreeCADUtils.createMesh(cs) for cs in cutSurfaces]

    pieces = FreeCADUtils.slicePart(body, cutSurfacesMesh)

    FreeCADUtils.addPartsToDocument(tabsSolid)
    FreeCADUtils.addPartsToDocument(leadsSolids)
    FreeCADUtils.addPartsToDocument(pieces)

    if os.path.exists("output/stl"):
        shutil.rmtree("output/stl", ignore_errors=True)

    os.makedirs("output/stl")

    for obj in document.Objects:
        filename = "output/stl/" + obj.Label + ".stl"
        obj.Shape.exportStl(filename)

    if not (checkParts(tabsSolid) and checkParts(leadsSolids) and checkParts(pieces)):
        raise SystemError

run()