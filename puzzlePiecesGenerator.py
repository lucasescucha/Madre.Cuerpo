import math
import os
import shutil
from curve.surfaceCurve import SurfaceCurve

import numpy as np
import pyclipper as clipper

from solid.operations import getLeadMesh, getShellMesh, getShellMeshFromSurfacesVertexs, surfaceToMeshSurface, surfaceToOffsetMeshSurface
from solver.offsetSolver import OffsetSolver
from surfaces.meshSurface import MeshSurface
from surfaces.offset.motherOffsetSurface import MotherOffsetSurface

import svg.utils as svgUtils
import surfaces.utils as sfcUtils
from sympy import sign
import utils.utils as utlUtils
import FreeCADTools.utils as FreeCADUtils

from configuration import configurationData
from surfaces.iSurface import ISurface

from surfaces.motherSurface import MotherSurface
from surfaces.shiftDecorator import ShiftDecorator
from surfaces.utils import calculateSamplingPoints, getUnitVector

from utils.utils import Configuration, rotateAroundAxis

TABS_HEIGHT_MARGIN = 0.1

X_VERSOR = np.array([1, 0, 0])
Y_VERSOR = np.array([0, 1, 0])
Z_VERSOR = np.array([0, 0, 1])

ZERO_VECTOR = ORIGIN = np.array([0, 0, 0])

SOLID_OPERATIONS_TOLERANCE = 0.1

def generateSurfaces(configuration: Configuration):
    def calculateReferenceSurface():
        surfaceConfig = configuration.surface
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

    thickness = configuration.manufacture.mold.puzzle.thickness
    tabsBaseThickness = configuration.manufacture.mold.puzzle.tabs.baseThickness

    baseOffset = thickness + tabsBaseThickness

    referenceSurface = calculateReferenceSurface()
    topOffsetSurface = MotherOffsetSurface(referenceSurface, thickness)
    baseOffsetSurface = MotherOffsetSurface(referenceSurface, baseOffset)

    return referenceSurface, topOffsetSurface, baseOffsetSurface

def generateSurfacesMeshes(configuration: Configuration, referenceSurface: ISurface):

    surface = configuration.surface.dimensions
    grid = configuration.manufacture.grid

    thickness = configuration.manufacture.mold.puzzle.thickness
    tabsBaseThickness = configuration.manufacture.mold.puzzle.tabs.baseThickness

    baseOffset = thickness + tabsBaseThickness

    xstart, xend = 0, surface.width/2
    xVect = calculateSamplingPoints(xstart, xend, referenceSurface, "x", grid.width)

    ystart, yend = (-surface.depth/2), (surface.depth / 2)
    yVect = calculateSamplingPoints(ystart, yend, referenceSurface, "y", grid.height)

    meshBottomSurface = surfaceToMeshSurface(referenceSurface, xVect, yVect)
    meshTopSurface = surfaceToOffsetMeshSurface(referenceSurface, thickness, xVect, yVect)
    meshBaseSurface = surfaceToOffsetMeshSurface(referenceSurface, baseOffset, xVect, yVect)

    surfaceMesh = getShellMesh(meshBottomSurface, meshTopSurface)
    baseMesh = getShellMesh(meshTopSurface, meshBaseSurface)

    return surfaceMesh, baseMesh

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
        polygon : np.array, direction, clearence, fastVersion = False):
    
    displacement = normal * dsplModule * (1 if direction == "up" else -1)

    if fastVersion:
        return polygon + 2*displacement

    comparer = lambda polp, sfcp: (polp > (sfcp+clearence) if direction == "up" else  polp < (sfcp-clearence))

    while any([not comparer(v[2], surface.F(v[0:2])) for v in polygon]):
        polygon = polygon + displacement
    
    return polygon

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

def generatePiecesCutSurfaces(configuration: Configuration, 
    referenceSurface : ISurface):

    def calculateTopAndBottomPoints(x, y, thickness):
        position2d = [x, y]
        position3d = np.append(position2d, [referenceSurface.F(position2d)]) 
        
        n_u = sfcUtils.getUnitNormalVector(referenceSurface, position2d)
        
        bottomOffset = getShiftedClearPoygon(referenceSurface, n_u, 3*thickness, 
                [position3d], "down", thickness)

        topOffset = getShiftedClearPoygon(referenceSurface, n_u, 3*thickness, 
                [position3d], "up", thickness)
            
        return [bottomOffset[0], topOffset[0]] 

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
        x = sfcUtils.arcLenght2Coordinate(referenceSurface, 
                    ix*pieceWidthLenght, "x", xstart)
        
        vertices = np.array(
                [calculateTopAndBottomPoints(x, y, thickness) for y in yVect])

        yield FreeCADUtils.createMesh(getLeadMesh(vertices[:, 0], vertices[:, 1]))

    for iy in range(1, yPieces):
        y = sfcUtils.arcLenght2Coordinate(referenceSurface, iy*pieceDepthLenght, "y", ystart)
        
        vertices = np.array(
                [calculateTopAndBottomPoints(x, y, thickness) for x in xVect])

        yield FreeCADUtils.createMesh(getLeadMesh(vertices[:, 0], vertices[:, 1]))

def generatePanelsLeads(configuration, referenceSurface):
    
    LEFT_LEAD, RIGHT_LEAD, TOP_LEAD, BOTTOM_LEAD = 0, 1, 2, 3

    TOP_FLANGE, BOTTOM_FLANGE, LEFT_FLANGE, RIGHT_FLANGE = 0, 1, 2, 3

    MIDDLE_FLANGE, CORNER_FLANGE = 0, 1

    TOP_DIRECTION, RIGHT_DIRECTION = 1, 1
    BOTTOM_DIRECTION, LEFT_DIRECTION = -1, -1

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
                t = getAxisTangentVector(p0, "y", TOP_DIRECTION)
            elif flangePosition == BOTTOM_FLANGE:
                t = getAxisTangentVector(p0, "x", LEFT_DIRECTION)
        
        if leadSide == RIGHT_LEAD:
            if flangePosition == TOP_FLANGE:
                t = getAxisTangentVector(p0, "x", RIGHT_DIRECTION)
            elif flangePosition == BOTTOM_FLANGE:
                t = getAxisTangentVector(p0, "y", BOTTOM_DIRECTION)
              
        if leadSide == TOP_LEAD:
            if flangePosition == LEFT_FLANGE:
                t = getAxisTangentVector(p0, "y", TOP_DIRECTION)
            elif flangePosition == RIGHT_FLANGE:
                t = getAxisTangentVector(p0, "x", RIGHT_DIRECTION)
        
        if leadSide == BOTTOM_LEAD:
            if flangePosition == LEFT_FLANGE:
                t = getAxisTangentVector(p0, "x", LEFT_DIRECTION)
            elif flangePosition == RIGHT_FLANGE:
                t = getAxisTangentVector(p0, "y", BOTTOM_DIRECTION)

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
        return np.array([referenceSurface.FOffset(p, r) for p in points])

    def getAxisTangentVector(p, axis, direction):
        grad = np.abs(referenceSurface.gradF(p))
        t = np.array([1, 0, grad[0]] if axis == "x" else [0, 1, grad[1]]) * direction
        t_hat = t/np.linalg.norm(t)
        
        return t_hat

    def extendVerticesWithVector(vertices, vector, distance):
        return np.array([(v + vector*distance) for v in vertices])

    def extendVertices(vertices, axis, direction, distance):
        return np.array([(v + (getAxisTangentVector(v[0:2], axis, direction)*distance)) 
                for v in vertices])

    def getSidePieceScrewsDrills(mainAxis, malstart, malend, salStart, 
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

            drillDirection = getAxisTangentVector(pos, secondaryAxis, saDir)
            drillCenter = getOffsetPoints([pos], moldThickness - (leadHeight/2))[0]

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

        v0 = getOffsetPoints([p0], moldThickness - (leadHeight/2))[0]

        d1, d2 = getFlangeDirections(leadSide, flangePosition, flangeType, p0)

        drillCenter = v0 + (d1 * flangeSize/2) - (d2 * leadThickness)
        drillDirection = d2 

        return FreeCADUtils.createCylinder(drillRadius, drillDepth, drillCenter, drillDirection)

    def getSideCoreMeshAndDrills(mainAxis, malStart, malEnd, salStart, saDir):
        def getSideCore(secondaryAxis, samplingPoints):
            mainAxisIndex = 0 if mainAxis == "x" else 1
            secAxisIndex = 1 if mainAxis == "x" else 0

            secAxisValue = samplingPoints[0, secAxisIndex]
            solver = OffsetSolver(
                SurfaceCurve(referenceSurface, mainAxis, secAxisValue))

            ar0Vertices = getOffsetPoints(samplingPoints, moldThickness)
            arefVertices = ar2Vertices = getOffsetPoints(samplingPoints, 0)
            
            ar3Vertices = np.empty_like(arefVertices)
            for i in range(len(samplingPoints)):
                sp = samplingPoints[i]
                offsetZ = solver.getOffsetFunctionZ(
                    [sp[mainAxisIndex]], -leadHeight)
                ar3Vertices[i, :] = np.append(sp, [offsetZ])

            aSideFace = [ar0Vertices, arefVertices, ar2Vertices, ar3Vertices]
            
            br0Vertices = extendVertices(ar0Vertices, secondaryAxis, saDir, leadThickness)
            brefVertices = extendVertices(arefVertices, secondaryAxis, saDir, leadThickness)
            
            vector = (X_VERSOR if secondaryAxis=="x" else Y_VERSOR) * saDir
            
            br2Vertices = extendVerticesWithVector(ar2Vertices, vector, leadThickness)
            br3Vertices = extendVerticesWithVector(ar3Vertices, vector, leadThickness)

            bSideFace = [br0Vertices, brefVertices, br2Vertices, br3Vertices]

            return getShellMeshFromSurfacesVertexs(aSideFace, bSideFace)

        secondaryAxis = "y" if mainAxis == "x" else "x"

        mainAxisIndex = 0 if mainAxis == "x" else 1
        secondaryAxisIndex = 0 if secondaryAxis == "x" else 1

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
        
        samplingPoints = np.zeros((len(mainAxisVect), 2))

        samplingPoints[:, mainAxisIndex] = mainAxisVect[:]
        samplingPoints[:, secondaryAxisIndex] = secondaryAxisPoint
        
        coreSideMesh = getSideCore(secondaryAxis, samplingPoints)
        sideDrills = getSidePieceScrewsDrills(mainAxis, malStart, malEnd, salStart, saDir, isPiece=False)

        return [[coreSideMesh, coreSideMesh], sideDrills]

    def getFlangeMeshAndDrill(leadSide, flangePosition, flangeType, lx, ly):
        if leadSide == LEFT_LEAD or leadSide == RIGHT_LEAD:
            if flangePosition == TOP_FLANGE:
                ly = ly + pieceDepthLenght
        else:
            if flangePosition == RIGHT_FLANGE:
                lx = lx + pieceWidthLenght

        p0 = sfcUtils.arcLenght2Coordinates(referenceSurface, [lx, ly], [xstart, ystart])

        rt = moldThickness
        rb = moldThickness - leadHeight

        d1, d2 = getFlangeDirections(leadSide, flangePosition, flangeType, p0)

        v0a = [getOffsetPoints([p0], rb)[0], getOffsetPoints([p0], rt)[0]]        
        v0b = v0a + d1 * flangeSize
        
        v1a = v0a + d2 * leadThickness
        v1b = v0b + d2 * leadThickness
        
        drill = getFlangeScrewDirll(leadSide, flangePosition, flangeType, lx, ly)

        return [getShellMeshFromSurfacesVertexs([v0a, v0b], [v1a, v1b]), [drill]]

    pnWidth = moldConfig.panels.dimensions.width
    pnHeight = moldConfig.panels.dimensions.height

    piecesDrills = []
    leadsAndDrills = []

    for ix in range(piecesX):        
        for iy in range(piecesY):
            rightLeadAndDrills, leftLeadAndDrills = [], []
            topLeadAndDrills, bottomLeadAndDrills = [], []

            leftSide = (ix % pnWidth) == 0
            rightSide = ((ix+1) % pnWidth) == 0

            bottomSide = (iy % pnHeight) == 0
            topSide = ((iy+1) % pnHeight) == 0
            
            lystart = iy*pieceDepthLenght
            lyend = lystart + pieceDepthLenght

            if leftSide:               
                lxstart = ix*pieceWidthLenght
                
                leftLeadAndDrills.append(
                    getSideCoreMeshAndDrills("y", lystart, lyend, lxstart, LEFT_DIRECTION))  

                piecesDrills.extend(
                    getSidePieceScrewsDrills("y", lystart, lyend, lxstart, LEFT_DIRECTION))     
            else:
                if topSide:
                    #agregar oreja a 90° arriba a la izquierda sólo a top 
                    topLeadAndDrills.append(
                        getFlangeMeshAndDrill(TOP_LEAD, LEFT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
                
                if bottomSide:
                    #agregar oreja a 90° abajo a la izquierda sólo a bottom 
                    bottomLeadAndDrills.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, LEFT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
            
            if rightSide:
                lxstart = (ix+1)*pieceWidthLenght

                rightLeadAndDrills.append(
                    getSideCoreMeshAndDrills("y", lystart, lyend, lxstart, RIGHT_DIRECTION))
                
                piecesDrills.extend(
                    getSidePieceScrewsDrills("y", lystart, lyend, lxstart, RIGHT_DIRECTION))     
            else:
                if topSide:
                    #agregar oreja a 90° arriba a la derecha sólo a top 
                    topLeadAndDrills.append(
                        getFlangeMeshAndDrill(TOP_LEAD, RIGHT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
                
                if bottomSide:
                    #agregar oreja a 90° abajo a la derecha sólo a bottom 
                    bottomLeadAndDrills.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, RIGHT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))

            lxstart = ix*pieceWidthLenght
            lxend = lxstart + pieceWidthLenght

            if topSide:
                lystart = (iy+1)*pieceDepthLenght

                topLeadAndDrills.append(
                    getSideCoreMeshAndDrills("x", lxstart, lxend, lystart, TOP_DIRECTION))
                
                piecesDrills.extend(
                    getSidePieceScrewsDrills("x", lxstart, lxend, lystart, TOP_DIRECTION))     
            else:
                if leftSide:
                    #agregar oreja a 90° arriba a la izquierda sólo a left
                    leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_LEAD, TOP_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 
                
                if rightSide:
                    #agregar oreja a 90° abajo a la derecha sólo a right 
                    rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_LEAD, BOTTOM_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 

            if bottomSide:
                lystart = iy*pieceDepthLenght

                bottomLeadAndDrills.append(
                    getSideCoreMeshAndDrills("x", lxstart, lxend, lystart, BOTTOM_DIRECTION))

                piecesDrills.extend(
                    getSidePieceScrewsDrills("x", lxstart, lxend, lystart, BOTTOM_DIRECTION))
            else:
                if leftSide:
                    #agregar oreja a 90° abajo a la izquierda sólo a left 
                    leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_LEAD, BOTTOM_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 
                
                if rightSide:
                    #agregar oreja a 90° arriba a la derecha sólo a right 
                    rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_LEAD, TOP_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 

            lxstart = ix*pieceWidthLenght
            lystart = iy*pieceDepthLenght

            if leftSide and topSide:
                #agregar oreja a 45° arriba a la izquierda a left y a top
                leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_LEAD, BOTTOM_FLANGE, CORNER_FLANGE, lxstart, lystart))
                topLeadAndDrills.append(
                        getFlangeMeshAndDrill(TOP_LEAD, LEFT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 
            
            if leftSide and bottomSide:
                #agregar oreja a 45° abajo a la izquierda a left y a bottom
                leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_LEAD, BOTTOM_FLANGE, CORNER_FLANGE, lxstart, lystart))
                bottomLeadAndDrills.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, LEFT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            if rightSide and topSide:
                #agregar oreja a 45° arriba a la derecha a right y a top
                rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_LEAD, TOP_FLANGE, CORNER_FLANGE, lxstart, lystart))
                topLeadAndDrills.append(
                        getFlangeMeshAndDrill(TOP_LEAD, RIGHT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            if rightSide and bottomSide:
                #agregar oreja a 45° abajo a la derecha a right y a bottom
                rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_LEAD, BOTTOM_FLANGE, CORNER_FLANGE, lxstart, lystart))
                bottomLeadAndDrills.append(
                        getFlangeMeshAndDrill(BOTTOM_LEAD, RIGHT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            def checkAndAppend(data):
                if len(data) > 0: leadsAndDrills.append(data)    
            
            checkAndAppend(rightLeadAndDrills)
            checkAndAppend(leftLeadAndDrills)
            checkAndAppend(topLeadAndDrills)
            checkAndAppend(bottomLeadAndDrills)

    return leadsAndDrills, piecesDrills

def checkParts(parts):
    for part in parts:
        if not part.isValid():
            return False

    return True

def run():
    configuration = Configuration(configurationData)

    if not configuration.check():
        raise ArithmeticError

    referenceSurface, topOffsetSurface, baseOffsetSurface = generateSurfaces(configuration)
    
    surfaceMesh, baseMesh = generateSurfacesMeshes(configuration, referenceSurface)

    document = FreeCADUtils.createNewDocument()

    surfaceSolid = FreeCADUtils.convertMeshToSolid(surfaceMesh)
    #baseSolid = FreeCADUtils.convertMeshToSolid(baseMesh)

    FreeCADUtils.addPartsToDocument([surfaceSolid])

    tabsShells = generateTabsShells(configuration, referenceSurface, topOffsetSurface, baseOffsetSurface)
    
    tabsShells = []
    tabsSolid = []

    body = surfaceSolid

    for shell in tabsShells:
        #FreeCADUtils.addMeshesToDocument([FreeCADUtils.createMesh(shell[0])])
        #result = FreeCADUtils.slicePart(body, [FreeCADUtils.createMesh(shell[0])])
        pass
        #volumeSortedResult = sorted(
        #    result, key=operator.attrgetter("Volume"), reverse=True)
        
        #body, tab = volumeSortedResult[0], volumeSortedResult[1:]

        #FreeCADUtils.addMeshesToDocument([FreeCADUtils.createMesh(shell[1])])
        #result = FreeCADUtils.slicePart(baseSolid, [FreeCADUtils.createMesh(shell[1])])

        #volumeSortedResult = sorted(
        #    result, key=operator.attrgetter("Volume"), reverse=True)

        #_, baseTab = volumeSortedResult[0], volumeSortedResult[1:]

        #tabsSolid.append(tab.fuse(baseTab))

    leadsAndDrills, piecesDrills = generatePanelsLeads(configuration, referenceSurface)
    for leadAndDrills in leadsAndDrills:
        for leadElement in leadAndDrills:
            if len(leadElement[0])==2:
                mesh = FreeCADUtils.convertMeshToSolid(leadElement[0][0])
                FreeCADUtils.addPartsToDocument([mesh])
                mesh = FreeCADUtils.convertMeshToSolid(leadElement[0][1])
                FreeCADUtils.addPartsToDocument([mesh])
            #FreeCADUtils.addPartsToDocument(leadElement[1])

    #for pieceDrills in piecesDrills:
        #body = body.cut(pieceDrills)
        #FreeCADUtils.addPartsToDocument([pieceDrills])

    surfaces = generatePiecesCutSurfaces(configuration, referenceSurface)
    
    FreeCADUtils.addMeshesToDocument(surfaces)

    

    #pieces = FreeCADUtils.slicePart(body, cutSurfacesMesh)

    #FreeCADUtils.addPartsToDocument(tabsSolid)
    #FreeCADUtils.addPartsToDocument(leadsSolids)
    #FreeCADUtils.addMeshesToDocument(cutSurfacesMesh)

    if not os.path.exists("output"):
        os.makedirs("output")
        
    if os.path.exists("output/output.FCStd"):
        os.remove("output/output.FCStd")

    FreeCADUtils.saveDocument(document, "output/output.FCStd")

    if os.path.exists("output/stl"):
        shutil.rmtree("output/stl", ignore_errors=True)

    #os.makedirs("output/stl")

    #for obj in document.Objects:
    #    filename = "output/stl/" + obj.Label + ".stl"
    #    obj.Shape.exportStl(filename)

    #if not (checkParts(tabsSolid) and checkParts(leadsSolids) and checkParts(pieces)):
    #    raise SystemError

run()