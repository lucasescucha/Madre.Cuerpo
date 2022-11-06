from generator.panelsWalls import generatePanelsWalls
from generator.pieces import generatePiecesCutSurfaces
from generator.surfaces import generateSurfaces, generateSurfacesMeshes
from generator.tabs import generateTabsShells
from utils.utils import Configuration

import FreeCADTools.utils as FreeCADUtils

def generateElements(configuration: Configuration):
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

    leadsAndDrills, piecesDrills = generatePanelsWalls(configuration, referenceSurface)
    for leadAndDrills in leadsAndDrills:
        for leadElement in leadAndDrills:
            mesh = FreeCADUtils.convertMeshToSolid(leadElement[0])
            FreeCADUtils.addPartsToDocument([mesh])

    #for pieceDrills in piecesDrills:
        #body = body.cut(pieceDrills)
        #FreeCADUtils.addPartsToDocument([pieceDrills])

    surfaces = generatePiecesCutSurfaces(configuration, referenceSurface)
    
    FreeCADUtils.addMeshesToDocument(surfaces)

    

    #pieces = FreeCADUtils.slicePart(body, cutSurfacesMesh)

    #FreeCADUtils.addPartsToDocument(tabsSolid)
    #FreeCADUtils.addPartsToDocument(leadsSolids)
    #FreeCADUtils.addMeshesToDocument(cutSurfacesMesh)

    return document


