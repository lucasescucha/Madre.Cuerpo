import numpy as np
import surfaces.utils as sfcUtils
import FreeCADTools.utils as FreeCADUtils

from generator.utils import calculatePiecesData
from utils.utils import Configuration, rotateAroundAxis
from solid.operations import getShellMeshFromSurfacesVertexs
from surfaces.iSurface import ISurface

def generatePanelsWalls(configuration: Configuration, referenceSurface: ISurface):
    
    LEFT_WALL, RIGHT_WALL, TOP_WALL, BOTTOM_WALL = 0, 1, 2, 3

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
        
        if leadSide == LEFT_WALL:
            if flangePosition == TOP_FLANGE:
                t = getAxisTangentVector(p0, "y", TOP_DIRECTION)
            elif flangePosition == BOTTOM_FLANGE:
                t = getAxisTangentVector(p0, "x", LEFT_DIRECTION)
        
        if leadSide == RIGHT_WALL:
            if flangePosition == TOP_FLANGE:
                t = getAxisTangentVector(p0, "x", RIGHT_DIRECTION)
            elif flangePosition == BOTTOM_FLANGE:
                t = getAxisTangentVector(p0, "y", BOTTOM_DIRECTION)
              
        if leadSide == TOP_WALL:
            if flangePosition == LEFT_FLANGE:
                t = getAxisTangentVector(p0, "y", TOP_DIRECTION)
            elif flangePosition == RIGHT_FLANGE:
                t = getAxisTangentVector(p0, "x", RIGHT_DIRECTION)
        
        if leadSide == BOTTOM_WALL:
            if flangePosition == LEFT_FLANGE:
                t = getAxisTangentVector(p0, "x", LEFT_DIRECTION)
            elif flangePosition == RIGHT_FLANGE:
                t = getAxisTangentVector(p0, "y", BOTTOM_DIRECTION)

        offsetNedeed = \
            (leadSide == LEFT_WALL and flangePosition == TOP_FLANGE) or \
            (leadSide == RIGHT_WALL and flangePosition == BOTTOM_FLANGE) or \
            (leadSide == TOP_WALL and flangePosition == RIGHT_FLANGE) or \
            (leadSide == BOTTOM_WALL and flangePosition == LEFT_FLANGE)
        
        if flangeType == MIDDLE_FLANGE:
            d1 = rotateAroundAxis(t, n_u, np.pi/2) if offsetNedeed else t
        elif flangeType == CORNER_FLANGE:
            d1 = rotateAroundAxis(t, n_u, np.pi/4)

        d2 = rotateAroundAxis(d1, n_u, np.pi/2 if offsetNedeed else -np.pi/2)

        return d1, d2

    def getOffsetPoints(points, r):
        return np.array([referenceSurface.FOffset(p, r) for p in points])

    def getAxisTangentVector(p, axis, direction):
        grad = referenceSurface.gradF(p)
        t = np.array([1, 0, grad[0]] if axis == "x" else [0, 1, grad[1]]) * direction
        
        return sfcUtils.getUnitVector(t)

    def extendVertices(vertices, axis, direction, distance):
        return np.array([(v + (getAxisTangentVector(v[0:2], axis, direction)*distance)) 
                for v in vertices])

    def getAxisParameters(leadSide, lx, ly):
        if (leadSide == TOP_WALL or leadSide == BOTTOM_WALL):
            mainAxis = "x"
            secondaryAxis = "y"
            malStart =  lx
            malEnd = malStart + pieceWidthLenght
            salStart = ly

            if leadSide == TOP_WALL:
                salStart += pieceDepthLenght
        else:
            mainAxis = "y"
            secondaryAxis = "x"
            malStart =  ly
            malEnd = malStart + pieceDepthLenght
            salStart = lx

            if leadSide == RIGHT_WALL:
                salStart += pieceWidthLenght
        
        return mainAxis, secondaryAxis, malStart, malEnd, salStart

    def getPieceWallScrewsDrills(leadSide, lx, ly, saDir, isPiece = True):

        mainAxis, secondaryAxis, malStart, malEnd, salStart = \
            getAxisParameters(leadSide, lx, ly)

        secondaryAxis = "y" if mainAxis == "x" else "x"

        pa_start = xstart if mainAxis == "x" else ystart
        sa_start = xstart if secondaryAxis == "x" else ystart

        sa_pos = sfcUtils.arcLenght2Coordinate(
                referenceSurface, salStart, secondaryAxis, sa_start)

        ldrills = (malEnd-malStart)/(leadScrews+1)
 
        for i in range(leadScrews):
            ldrill = malStart + (i+1)*ldrills 
            
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

    def getCoreWallMeshAndDrills(leadSide, lx, ly, saDir):
        
        mainAxis, secondaryAxis, malStart, malEnd, salStart = \
            getAxisParameters(leadSide, lx, ly)

        pa_start = xstart if mainAxis == "x" else ystart
        sa_start = xstart if secondaryAxis == "x" else ystart
        
        gridDim = grid.width if mainAxis == "x" else grid.height

        mpa_start = sfcUtils.arcLenght2Coordinate(
                referenceSurface, malStart, mainAxis, pa_start)
        mpa_end = sfcUtils.arcLenght2Coordinate(
                referenceSurface, malEnd, mainAxis, pa_start)
        
        secondaryAxisPoint = sfcUtils.arcLenght2Coordinate(
                referenceSurface, salStart, secondaryAxis, sa_start)

        mainAxisVect = sfcUtils.calculateSamplingPoints(mpa_start, mpa_end, 
                            referenceSurface, mainAxis, gridDim)

        if mainAxis == "x":
            points = [[x, secondaryAxisPoint] for x in mainAxisVect]
        else:
            points = [[secondaryAxisPoint, y] for y in mainAxisVect]

        rt = moldThickness
        rb = moldThickness - leadHeight

        topVertices = getOffsetPoints(points, rt)
        bottomVertices = getOffsetPoints(points, rb)

        leadAFace = [topVertices, bottomVertices]
        
        eTopVertices = extendVertices(topVertices, secondaryAxis, saDir, leadThickness)
        eBottomVertices = extendVertices(bottomVertices, secondaryAxis, saDir, leadThickness)

        leadBFace = [eTopVertices, eBottomVertices]

        drills = getPieceWallScrewsDrills(leadSide, lx, ly, saDir, isPiece=False)

        return [getShellMeshFromSurfacesVertexs(leadAFace, leadBFace), drills]

    def getFlangeMeshAndDrill(leadSide, flangePosition, flangeType, lx, ly):
        if leadSide == TOP_WALL or flangePosition == TOP_FLANGE:
            ly = ly + pieceDepthLenght
        
        if leadSide == RIGHT_WALL or flangePosition == RIGHT_FLANGE:
            lx = lx + pieceWidthLenght
        
        axisIndex = 1 if (leadSide == LEFT_WALL or leadSide == RIGHT_WALL) else 0

        p0 = sfcUtils.arcLenght2Coordinates(referenceSurface, [lx, ly], [xstart, ystart])

        rt = moldThickness
        rb = moldThickness - leadHeight

        d1, d2 = getFlangeDirections(leadSide, flangePosition, flangeType, p0)

        v0a = [getOffsetPoints([p0], rb)[0], getOffsetPoints([p0], rt)[0]]        
        v0b = v0a + d1 * flangeSize
        
        if flangeType == MIDDLE_FLANGE:
            p0[axisIndex] = p0[axisIndex] + (d2 * leadThickness)[axisIndex]
            v1a = [getOffsetPoints([p0], rb)[0], getOffsetPoints([p0], rt)[0]]        
        else:
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
            lxstart = ix*pieceWidthLenght

            if leftSide:               
                leftLeadAndDrills.append(
                    getCoreWallMeshAndDrills(LEFT_WALL, lxstart, lystart, LEFT_DIRECTION))  

                piecesDrills.extend(
                    getPieceWallScrewsDrills(LEFT_WALL, lxstart, lystart, LEFT_DIRECTION))     
            else:
                if topSide:
                    #agregar oreja a 90° arriba a la izquierda sólo a top 
                    topLeadAndDrills.append(
                        getFlangeMeshAndDrill(TOP_WALL, LEFT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
                
                if bottomSide:
                    #agregar oreja a 90° abajo a la izquierda sólo a bottom 
                    bottomLeadAndDrills.append(
                        getFlangeMeshAndDrill(BOTTOM_WALL, LEFT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
            
            if rightSide:
                rightLeadAndDrills.append(
                    getCoreWallMeshAndDrills(RIGHT_WALL, lxstart, lystart, RIGHT_DIRECTION))
                
                piecesDrills.extend(
                    getPieceWallScrewsDrills(RIGHT_WALL, lxstart, lystart, RIGHT_DIRECTION))     
            else:
                if topSide:
                    #agregar oreja a 90° arriba a la derecha sólo a top 
                    topLeadAndDrills.append(
                        getFlangeMeshAndDrill(TOP_WALL, RIGHT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))
                
                if bottomSide:
                    #agregar oreja a 90° abajo a la derecha sólo a bottom 
                    bottomLeadAndDrills.append(
                        getFlangeMeshAndDrill(BOTTOM_WALL, RIGHT_FLANGE, MIDDLE_FLANGE, lxstart, lystart))

            if topSide:
                topLeadAndDrills.append(
                    getCoreWallMeshAndDrills(TOP_WALL, lxstart, lystart, TOP_DIRECTION))
                
                piecesDrills.extend(
                    getPieceWallScrewsDrills(TOP_WALL, lxstart, lystart, TOP_DIRECTION))     
            else:
                if leftSide:
                    #agregar oreja a 90° arriba a la izquierda sólo a left
                    leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_WALL, TOP_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 
                
                if rightSide:
                    #agregar oreja a 90° abajo a la derecha sólo a right 
                    rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_WALL, BOTTOM_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 

            if bottomSide:
                bottomLeadAndDrills.append(
                    getCoreWallMeshAndDrills(BOTTOM_WALL, lxstart, lystart, BOTTOM_DIRECTION))

                piecesDrills.extend(
                    getPieceWallScrewsDrills(BOTTOM_WALL, lxstart, lystart, BOTTOM_DIRECTION))
            else:
                if leftSide:
                    #agregar oreja a 90° abajo a la izquierda sólo a left 
                    leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_WALL, BOTTOM_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 
                
                if rightSide:
                    #agregar oreja a 90° arriba a la derecha sólo a right 
                    rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_WALL, TOP_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 

            if leftSide and topSide:
                #agregar oreja a 45° arriba a la izquierda a left y a top
                leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_WALL, TOP_FLANGE, CORNER_FLANGE, lxstart, lystart))
                topLeadAndDrills.append(
                        getFlangeMeshAndDrill(TOP_WALL, LEFT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 
            
            if leftSide and bottomSide:
                #agregar oreja a 45° abajo a la izquierda a left y a bottom
                leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_WALL, BOTTOM_FLANGE, CORNER_FLANGE, lxstart, lystart))
                bottomLeadAndDrills.append(
                        getFlangeMeshAndDrill(BOTTOM_WALL, LEFT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            if rightSide and topSide:
                #agregar oreja a 45° arriba a la derecha a right y a top
                rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_WALL, TOP_FLANGE, CORNER_FLANGE, lxstart, lystart))
                topLeadAndDrills.append(
                        getFlangeMeshAndDrill(TOP_WALL, RIGHT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            if rightSide and bottomSide:
                #agregar oreja a 45° abajo a la derecha a right y a bottom
                rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_WALL, BOTTOM_FLANGE, CORNER_FLANGE, lxstart, lystart))
                bottomLeadAndDrills.append(
                        getFlangeMeshAndDrill(BOTTOM_WALL, RIGHT_FLANGE, CORNER_FLANGE, lxstart, lystart)) 

            def checkAndAppend(data):
                if len(data) > 0: leadsAndDrills.append(data)    
            
            checkAndAppend(rightLeadAndDrills)
            checkAndAppend(leftLeadAndDrills)
            checkAndAppend(topLeadAndDrills)
            checkAndAppend(bottomLeadAndDrills)

    return leadsAndDrills, piecesDrills