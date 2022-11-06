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

    def getFlangeDirections(wallSide, flangePosition, flangeType, p):
        n_u = sfcUtils.getUnitNormalVector(referenceSurface, p)
        
        if flangeType == CORNER_FLANGE:
            if wallSide == LEFT_WALL:
                if flangePosition == TOP_FLANGE:
                    t = getAxisTangentVector(p, "y", TOP_DIRECTION)
                elif flangePosition == BOTTOM_FLANGE:
                    t = getAxisTangentVector(p, "x", LEFT_DIRECTION)
            
            #
            if wallSide == RIGHT_WALL:
                if flangePosition == TOP_FLANGE:
                    ###
                    e1 = getExtensionVector(TOP_WALL, p)
                    e2 = getExtensionVector(RIGHT_WALL, p)
                    t = sfcUtils.getUnitVector(e1+e2)
                    #t = getAxisTangentVector(p, "x", RIGHT_DIRECTION)
                elif flangePosition == BOTTOM_FLANGE:
                    t = getAxisTangentVector(p, "y", BOTTOM_DIRECTION)
                
            if wallSide == TOP_WALL:
                if flangePosition == LEFT_FLANGE:
                    t = getAxisTangentVector(p, "y", TOP_DIRECTION)
                elif flangePosition == RIGHT_FLANGE:
                    ####
                    d1 = getExtensionVector(TOP_WALL, p)
                    d1 = getExtensionVector(RIGHT_WALL, p)
                    #t = getAxisTangentVector(p, "x", RIGHT_DIRECTION)
            
            if wallSide == BOTTOM_WALL:
                if flangePosition == LEFT_FLANGE:
                    t = getAxisTangentVector(p, "x", LEFT_DIRECTION)
                elif flangePosition == RIGHT_FLANGE:
                    t = getAxisTangentVector(p, "y", BOTTOM_DIRECTION)

            offsetNedeed = \
                (wallSide == LEFT_WALL and flangePosition == TOP_FLANGE) or \
                (wallSide == RIGHT_WALL and flangePosition == BOTTOM_FLANGE) or \
                (wallSide == TOP_WALL and flangePosition == RIGHT_FLANGE) or \
                (wallSide == BOTTOM_WALL and flangePosition == LEFT_FLANGE)
            
            d1 = rotateAroundAxis(t, n_u, np.pi/4)
            d2 = rotateAroundAxis(d1, n_u, np.pi/2 if offsetNedeed else -np.pi/2)

            return d1, d2
        else:
            d1 = getExtensionVector(wallSide, p)

            if wallSide == LEFT_WALL:
                if flangePosition == TOP_FLANGE:
                    d2 = rotateAroundAxis(d1, n_u, np.pi/2)
                elif flangePosition == BOTTOM_FLANGE:
                    d2 = rotateAroundAxis(d1, n_u, -np.pi/2)
            
            if wallSide == RIGHT_WALL:
                if flangePosition == TOP_FLANGE:
                    d2 = rotateAroundAxis(d1, n_u, -np.pi/2)
                elif flangePosition == BOTTOM_FLANGE:
                    d2 = rotateAroundAxis(d1, n_u, np.pi/2)
                
            if wallSide == TOP_WALL:
                if flangePosition == LEFT_FLANGE:
                    d2 = rotateAroundAxis(d1, n_u, -np.pi/2)
                elif flangePosition == RIGHT_FLANGE:
                    d2 = rotateAroundAxis(d1, n_u, np.pi/2)
            
            if wallSide == BOTTOM_WALL:
                if flangePosition == LEFT_FLANGE:
                    d2 = rotateAroundAxis(d1, n_u, np.pi/2)
                elif flangePosition == RIGHT_FLANGE:
                    d2 = rotateAroundAxis(d1, n_u, -np.pi/2)

            return d1, d2

    def getOffsetPoints(points, r):
        return np.array([referenceSurface.FOffset(p, r) for p in points])

    def getAxisTangentVector(p, axis, direction):
        grad = referenceSurface.gradF(p)
        t = np.array([1, 0, grad[0]] if axis == "x" else [0, 1, grad[1]]) * direction
        
        return sfcUtils.getUnitVector(t)

    def getExtensionVector(wallSide, p):
        n_u = sfcUtils.getUnitNormalVector(referenceSurface, p)
        if wallSide == BOTTOM_WALL:
            t = getAxisTangentVector(p, "x", LEFT_DIRECTION)
        elif wallSide == TOP_WALL:
            t = getAxisTangentVector(p, "x", RIGHT_DIRECTION)
        elif wallSide == LEFT_WALL:
            t = getAxisTangentVector(p, "y", TOP_DIRECTION)
        elif wallSide == RIGHT_WALL:
            t = getAxisTangentVector(p, "y", BOTTOM_DIRECTION)
            
        return rotateAroundAxis(t, n_u, np.pi/2)
        
    def getAxisParameters(wallSide, lx, ly):
        if (wallSide == TOP_WALL or wallSide == BOTTOM_WALL):
            mainAxis = "x"
            secondaryAxis = "y"
            malStart =  lx
            malEnd = malStart + pieceWidthLenght
            salStart = ly

            if wallSide == TOP_WALL:
                salStart += pieceDepthLenght
        else:
            mainAxis = "y"
            secondaryAxis = "x"
            malStart =  ly
            malEnd = malStart + pieceDepthLenght
            salStart = lx

            if wallSide == RIGHT_WALL:
                salStart += pieceWidthLenght
        
        return mainAxis, secondaryAxis, malStart, malEnd, salStart

    def getPieceWallScrewsDrills(wallSide, lx, ly, isPiece = True):

        mainAxis, secondaryAxis, malStart, malEnd, salStart = \
            getAxisParameters(wallSide, lx, ly)

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

            drillDirection = getExtensionVector(wallSide, pos)
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

    def getFlangeScrewDirll(wallSide, flangePosition, flangeType, lx, ly):
        drillDepth = 3*leadThickness
        drillRadius = leadScrewsDiameter/2

        p0 = sfcUtils.arcLenght2Coordinates(referenceSurface, [lx, ly], [xstart, ystart])

        v0 = getOffsetPoints([p0], moldThickness - (leadHeight/2))[0]

        d1, d2 = getFlangeDirections(wallSide, flangePosition, flangeType, p0)

        drillCenter = v0 + (d1 * flangeSize/2) - (d2 * leadThickness)
        drillDirection = d2 

        return FreeCADUtils.createCylinder(drillRadius, drillDepth, drillCenter, drillDirection)

    def getCoreWallMeshAndDrills(wallSide, lx, ly):
        def extendVertices(samplingPoint, vertices):
            for i in range(len(samplingPoint)):
                yield vertices[i] + getExtensionVector(wallSide, samplingPoint[i])*leadThickness

        mainAxis, secondaryAxis, malStart, malEnd, salStart = \
            getAxisParameters(wallSide, lx, ly)

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
            samplingPoints = list([[x, secondaryAxisPoint] for x in mainAxisVect])
        else:
            samplingPoints = list([[secondaryAxisPoint, y] for y in mainAxisVect])

        rt = moldThickness
        rb = moldThickness - leadHeight

        topVertices = getOffsetPoints(samplingPoints, rt)
        bottomVertices = getOffsetPoints(samplingPoints, rb)

        leadAFace = [topVertices, bottomVertices]
        
        eTopVertices = list(extendVertices(samplingPoints, topVertices))
        eBottomVertices = list(extendVertices(samplingPoints, bottomVertices))

        leadBFace = [eTopVertices, eBottomVertices]

        drills = getPieceWallScrewsDrills(wallSide, lx, ly, isPiece=False)

        return [getShellMeshFromSurfacesVertexs(leadAFace, leadBFace), drills]

    def getFlangeMeshAndDrill(wallSide, flangePosition, flangeType, lx, ly):
        if wallSide == TOP_WALL or flangePosition == TOP_FLANGE:
            ly = ly + pieceDepthLenght
        
        if wallSide == RIGHT_WALL or flangePosition == RIGHT_FLANGE:
            lx = lx + pieceWidthLenght
        
        axisIndex = 1 if (wallSide == LEFT_WALL or wallSide == RIGHT_WALL) else 0

        p0 = sfcUtils.arcLenght2Coordinates(referenceSurface, [lx, ly], [xstart, ystart])

        rt = moldThickness
        rb = moldThickness - leadHeight

        d1, d2 = getFlangeDirections(wallSide, flangePosition, flangeType, p0)

        v0a = [getOffsetPoints([p0], rb)[0], getOffsetPoints([p0], rt)[0]]        
        v0b = v0a + d1 * flangeSize
        
        if flangeType == MIDDLE_FLANGE:
            p0[axisIndex] = p0[axisIndex] + (d2 * leadThickness)[axisIndex]
            v1a = [getOffsetPoints([p0], rb)[0], getOffsetPoints([p0], rt)[0]]        
        else:
            v1a = v0a + d2 * leadThickness

        v1b = v0b + d2 * leadThickness
        
        drill = getFlangeScrewDirll(wallSide, flangePosition, flangeType, lx, ly)

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
                    getCoreWallMeshAndDrills(LEFT_WALL, lxstart, lystart))  

                piecesDrills.extend(
                    getPieceWallScrewsDrills(LEFT_WALL, lxstart, lystart))     
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
                    getCoreWallMeshAndDrills(RIGHT_WALL, lxstart, lystart))
                
                piecesDrills.extend(
                    getPieceWallScrewsDrills(RIGHT_WALL, lxstart, lystart))     
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
                    getCoreWallMeshAndDrills(TOP_WALL, lxstart, lystart))
                
                piecesDrills.extend(
                    getPieceWallScrewsDrills(TOP_WALL, lxstart, lystart))     
            else:
                if leftSide:
                    #agregar oreja a 90° arriba a la izquierda sólo a left
                    leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_WALL, TOP_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 
                
                if rightSide:
                    #agregar oreja a 90° abajo a la derecha sólo a right 
                    rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_WALL, TOP_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 

            if bottomSide:
                bottomLeadAndDrills.append(
                    getCoreWallMeshAndDrills(BOTTOM_WALL, lxstart, lystart))

                piecesDrills.extend(
                    getPieceWallScrewsDrills(BOTTOM_WALL, lxstart, lystart))
            else:
                if leftSide:
                    #agregar oreja a 90° abajo a la izquierda sólo a left 
                    leftLeadAndDrills.append(
                        getFlangeMeshAndDrill(LEFT_WALL, BOTTOM_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 
                
                if rightSide:
                    #agregar oreja a 90° arriba a la derecha sólo a right 
                    rightLeadAndDrills.append(
                        getFlangeMeshAndDrill(RIGHT_WALL, BOTTOM_FLANGE, MIDDLE_FLANGE, lxstart, lystart)) 

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