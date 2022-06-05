import operator
import os
import shutil

import config

import numpy as np

import FreeCADTools.utils as FreeCADUtils
import svg.utils as svgUtils

from surfaces.iSurface import ISurface

from surfaces.offset.motherOffsetSurface import MotherOffsetSurface
from surfaces.motherSurface import MotherSurface
from surfaces.shiftDecorator import ShiftDecorator

from solid.operations import getShellMesh

from utils.utils import Configuration

TABS_HEIGHT_MARGIN = 0.1


def checkConfguration(configuration: Configuration):
    def checkRatio(a, b):
        return (a/b).is_integer()

    sDimensionsConfig = configuration.surface.dimensions
    gridConfig = configuration.manufacture.grid
    moldConfig = configuration.manufacture.mold
    puzzleConfig = moldConfig.puzzle

    result = True

    result &= checkRatio(sDimensionsConfig.width/2, gridConfig.width)
    result &= checkRatio(sDimensionsConfig.depth, gridConfig.height)

    result &= checkRatio(sDimensionsConfig.width/2, puzzleConfig.pieces.width)
    result &= checkRatio(sDimensionsConfig.depth, puzzleConfig.pieces.height)

    piecesX = int((sDimensionsConfig.width/2)/puzzleConfig.pieces.width)
    piecesY = int(sDimensionsConfig.depth/puzzleConfig.pieces.height)

    result &= checkRatio(piecesX, moldConfig.panels.dimensions.width)
    result &= checkRatio(piecesY, moldConfig.panels.dimensions.height)

    return result

def generateTopBottomSurfaces(configuration: Configuration):
    def calculateSurfaceParmeters(surfaceConfig):
        zOffset = surfaceConfig.shape.zeroHeight

        # z = zOffset + a*x^2 -  b*y^2

        # y=0 -> z = zOffset + a*x^2 -> height = zOffset + a*(width/2)^2
        # a = (height - zOffset) / (width/2)^2
        a = (surfaceConfig.dimensions.height - zOffset) / \
            ((surfaceConfig.dimensions.width/2)**2)

        # x=0 -> z = zOffset - b*y^2 -> 0 = zOffset - b*(depth/2)^2
        # b = zOffset / (depth/2)^2
        b = zOffset / ((surfaceConfig.dimensions.depth/2)**2)

        return a, b, zOffset

    a, b, zOffset = calculateSurfaceParmeters(configuration.surface)

    bottomSurface = ShiftDecorator(
        MotherSurface(a, b),  np.array([0, 0, zOffset]))

    thickness = configuration.manufacture.mold.puzzle.thickness

    topSurface = MotherOffsetSurface(bottomSurface, thickness)

    return topSurface, bottomSurface

def generateSurfaceSolid(configuration: Configuration,
                         topSurface: ISurface, bottomSurface: ISurface):
    dimensions = configuration.surface.dimensions
    grid = configuration.manufacture.grid

    # Because of the surface's symmetry, xstart=0
    xstart, xend, xsteps = 0, dimensions.width/2, \
        int((dimensions.width/2)/grid.width)

    ystart, yend, ysteps = (-dimensions.depth/2), (dimensions.depth / 2), \
        int(dimensions.depth/grid.height)

    shellTrianglesMesh = getShellMesh(topSurface, bottomSurface,
                                      xstart, xend, xsteps, ystart, yend, ysteps)

    return FreeCADUtils.convertMeshToSolid(shellTrianglesMesh)

def generateTabsShells(configuration: Configuration,
                       topSurface: ISurface, bottomSurface: ISurface):
    
    def getExtrudedPolygon(topSurface, bottomSurface, tabPolygon):
        def getMaxMinValuesShell(topSurface, bottomSurface, points):
            def evalSurface(surface, points, f):
                return f([surface.F(p) for p in points])

            min = evalSurface(bottomSurface, points, np.min)
            max = evalSurface(topSurface, points,  np.max)

            return min, max

        min, max = getMaxMinValuesShell(topSurface, bottomSurface, tabPolygon)

        zStart = min-(max-min)*TABS_HEIGHT_MARGIN
        zEnd = max + (max-min)*TABS_HEIGHT_MARGIN

        return FreeCADUtils.extrudePolygonInZ(tabPolygon, zEnd, zStart)

    surfDims = configuration.surface.dimensions
    moldConfig = configuration.manufacture.mold
    puzzleConfig = moldConfig.puzzle
    pDimensionsConfig = moldConfig.panels.dimensions

    sWidth, sDepth = surfDims.width, surfDims.depth
    pWidth, pHeight = puzzleConfig.pieces.width, puzzleConfig.pieces.height
    pnlWidth, pnlHeight = pDimensionsConfig.width, pDimensionsConfig.height

    puzzleTabPolygon = svgUtils.getTabJoinPolygonFromSVGFile(
        puzzleConfig.tabs.file, puzzleConfig.tabs.width)

    piecesX = int((sWidth/2)/pWidth)
    piecesY = int(sDepth/pHeight)

    for ix in range(int(piecesX)):
        for iy in range(1, int(piecesY)):

            if iy % pnlHeight == 0:
                continue

            position = np.array([(ix+1/2)*pWidth, iy*pHeight - (sDepth/2)])

            tabPolygon = puzzleTabPolygon + position

            yield getExtrudedPolygon(topSurface, bottomSurface, tabPolygon)

    R = np.array(((0, -1), (1,  0)))
    rotPuzzleTabPolygon = np.dot(puzzleTabPolygon, R.T)

    for iy in range(int(piecesY)):
        for ix in range(1, int(piecesX)):

            if ix % pnlWidth == 0:
                continue

            position = np.array([ix*pWidth, (iy+1/2)*pHeight - (sDepth/2)])

            tabPolygon = rotPuzzleTabPolygon + position

            yield getExtrudedPolygon(topSurface, bottomSurface, tabPolygon)

def generatePiecesCutPlanes(configuration: Configuration):
    surfDims = configuration.surface.dimensions
    puzzleConfig = configuration.manufacture.mold.puzzle

    pWidth, pHeight = puzzleConfig.pieces.width, puzzleConfig.pieces.height
    sWidth, sDepth, sHight = surfDims.width, surfDims.depth, surfDims.height*2  # FIX IT!!

    piecesX = int((sWidth/2)/pWidth)
    piecesY = int(sDepth/pHeight)

    for ix in range(1, piecesX):
        pos = [ix*pWidth, sDepth/2, 0]
        yield FreeCADUtils.createPlane(sHight, sDepth, pos, [1, 0, 0])

    for iy in range(1, piecesY):
        pos = [0, iy*pHeight - (sDepth/2), 0]
        yield FreeCADUtils.createPlane(sHight, sWidth/2, pos, [0, 1, 0])

def generateParts(configuration, topSurface, bottomSurface):
    surfaceSolid = generateSurfaceSolid(
        configuration, topSurface, bottomSurface)

    tabsShells = generateTabsShells(
        configuration, topSurface, bottomSurface)

    return surfaceSolid, list(tabsShells)

def generatePanelsLeads(configuration, topSurface, bottomSurface):
    surfDims = configuration.surface.dimensions
    moldConfig = configuration.manufacture.mold
    puzzleConfig = moldConfig.puzzle

    moldThickness = moldConfig.puzzle.thickness

    leadHeight = moldConfig.panels.leads.height
    leadThickness = moldConfig.panels.leads.thickness

    leadBottomSurface = MotherOffsetSurface(
        bottomSurface, moldThickness - leadHeight)

    grid = configuration.manufacture.grid

    pWidth, pHeight = puzzleConfig.pieces.width, puzzleConfig.pieces.height
    sWidth, sDepth = surfDims.width, surfDims.depth

    piecesX = int((sWidth/2)/pWidth)
    piecesY = int(sDepth/pHeight)

    pnWidth = moldConfig.panels.dimensions.width
    pnHeight = moldConfig.panels.dimensions.height

    def getLeadPart(xstart, xend, xsteps, ystart, yend, ysteps):
        leadTriMesh = getShellMesh(topSurface, leadBottomSurface,
                                   xstart, xend, xsteps, ystart, yend, ysteps)
        return FreeCADUtils.convertMeshToSolid(leadTriMesh)

    leads = []
    for ix in range(piecesX):
        pieceLeads = []
        for iy in range(piecesY):
            faceLeads = []

            leftSide = (ix % pnWidth) == 0
            rightSide = ((ix+1) % pnWidth) == 0

            topSide = (iy % pnHeight) == 0
            bottomSide = ((iy+1) % pnHeight) == 0

            # Side leads
            xsteps = 1

            ystart, yend = (iy*pHeight)-(sDepth/2), ((iy+1)*pHeight)-(sDepth/2)
            ysteps = int(round((yend - ystart)/grid.height))

            if leftSide:
                xend = ix*pWidth
                xstart = xend - leadThickness

                faceLeads.append(getLeadPart(
                    xstart, xend, xsteps, ystart, yend, ysteps))

            if rightSide:
                xstart = (ix+1)*pWidth
                xend = xstart + leadThickness

                faceLeads.append(getLeadPart(
                    xstart, xend, xsteps, ystart, yend, ysteps))

            ysteps = 1

            xstart, xend = ix*pWidth, (ix+1)*pWidth
            xsteps = int(round((xend - xstart)/grid.width))

            if topSide:
                yend = (iy*pHeight)-(sDepth/2)
                ystart = yend - leadThickness

                faceLeads.append(getLeadPart(
                    xstart, xend, xsteps, ystart, yend, ysteps))

            if bottomSide:
                ystart = ((iy+1)*pHeight)-(sDepth/2)
                yend = ystart + leadThickness

                faceLeads.append(getLeadPart(
                    xstart, xend, xsteps, ystart, yend, ysteps))

            # Corners
            xsteps, ysteps = 1, 1

            xend = ix*pWidth
            xstart = xend - leadThickness

            if leftSide and topSide:
                yend = (iy*pHeight)-(sDepth/2)
                ystart = yend - leadThickness

                faceLeads.append(getLeadPart(
                    xstart, xend, xsteps, ystart, yend, ysteps))

            if leftSide and bottomSide:
                ystart = ((iy+1)*pHeight)-(sDepth/2)
                yend = ystart + leadThickness

                faceLeads.append(getLeadPart(
                    xstart, xend, xsteps, ystart, yend, ysteps))

            xstart = (ix+1)*pWidth
            xend = xstart + leadThickness

            if rightSide and topSide:
                yend = (iy*pHeight)-(sDepth/2)
                ystart = yend - leadThickness

                faceLeads.append(getLeadPart(
                    xstart, xend, xsteps, ystart, yend, ysteps))

            if rightSide and bottomSide:
                ystart = ((iy+1)*pHeight)-(sDepth/2)
                yend = ystart + leadThickness

                faceLeads.append(getLeadPart(
                    xstart, xend, xsteps, ystart, yend, ysteps))

            pieceLeads.append(faceLeads)
        leads.append(pieceLeads)

    return leads

def sliceTabsFromSurface(surfaceSolid, tabsShells):
    result = FreeCADUtils.slicePart(surfaceSolid, tabsShells)
    volumeSortedResult = sorted(
        result, key=operator.attrgetter("Volume"), reverse=True)

    return volumeSortedResult[0], volumeSortedResult[1:]

def sliceTabsAndPieces(surfaceSolid, tabsShells, planes):
    body, tabs = sliceTabsFromSurface(surfaceSolid, tabsShells)
    pieces = FreeCADUtils.slicePart(body, planes)

    return tabs, pieces

def joinPiecesLeads(configuration, orderedPieces, leads):
    surfDims = configuration.surface.dimensions
    puzzleConfig = configuration.manufacture.mold.puzzle

    pWidth, pHeight = puzzleConfig.pieces.width, puzzleConfig.pieces.height
    sWidth, sDepth = surfDims.width, surfDims.depth

    piecesX = int((sWidth/2)/pWidth)
    piecesY = int(sDepth/pHeight)

    parts = []
    for ix in range(piecesX):
        for iy in range(piecesY):
            part = orderedPieces[ix][iy]
            # parts.append(part)
            for lead in leads[ix][iy]:
                part = part.fuse(lead)
                # parts.append(lead)

            parts.append(part)

    return parts

def checkParts(parts):
    for part in parts:
        if not part.isValid():
            return False

    return True

def orderElementInGrid(gridWidth, gridHeight, elements, fGridPoint, fElementPosition):
    def dist(elementPos, gridPos): return np.linalg.norm(elementPos - gridPos)

    orderedElements = []
    for ix in range(gridWidth):
        rowElements = []
        for iy in range(gridHeight):

            distances = [
                [e, dist(fElementPosition(e), fGridPoint(ix, iy))] for e in elements]

            nearestPiece = (min(distances, key=lambda e: e[1]))[0]
            rowElements.append(nearestPiece)

        orderedElements.append(rowElements)
    return orderedElements

def orderPieces(configuration, pieces, bottomSurface):
    surfDims = configuration.surface.dimensions
    puzzleConfig = configuration.manufacture.mold.puzzle

    pWidth, pHeight = puzzleConfig.pieces.width, puzzleConfig.pieces.height
    sWidth, sDepth = surfDims.width, surfDims.depth

    piecesX = int((sWidth/2)/pWidth)
    piecesY = int(sDepth/pHeight)

    def fPiecePosition(piece): return np.array(piece.CenterOfMass)

    def fGridPoint(ix, iy):
        x = (ix+1/2)*pWidth
        y = (iy+1/2)*pHeight - (sDepth/2)
        z = bottomSurface.F(np.array([x, y]))

        return np.array([x, y, z])

    return orderElementInGrid(piecesX, piecesY, pieces, fGridPoint, fPiecePosition)

def run():
    configuration = Configuration(config.configuration)

    if not checkConfguration(configuration):
        raise ArithmeticError

    topSurface, bottomSurface = generateTopBottomSurfaces(configuration)

    document = FreeCADUtils.createNewDocument()

    surfaceSolid, tabsShells = generateParts(
        configuration, topSurface, bottomSurface)

    planes = list(generatePiecesCutPlanes(configuration))

    tabs, pieces = sliceTabsAndPieces(surfaceSolid, tabsShells, planes)

    orderedPieces = orderPieces(configuration, pieces, bottomSurface)

    leads = generatePanelsLeads(configuration, topSurface, bottomSurface)

    pieces = joinPiecesLeads(configuration, orderedPieces, leads)

    FreeCADUtils.addPartsToDocument(tabs)
    FreeCADUtils.addPartsToDocument(pieces)

    if not os.path.exists("output"):
        os.makedirs("output")
        
    if os.path.exists("output/output.FCStd"):
        os.remove("output/output.FCStd")

    FreeCADUtils.saveDocument(document, "output/output.FCStd")

    if os.path.exists("output/stl"):
        shutil.rmtree("output/stl", ignore_errors=True)

    os.makedirs("output/stl")

    for obj in document.Objects:
        filename = "output/stl/" + obj.Label + ".stl"
        obj.Shape.exportStl(filename)

    if not (checkParts(tabs) and checkParts(pieces)):
        raise SystemError

run()