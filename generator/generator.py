from decimal import InvalidOperation
from generator.panelsWalls import generatePanelsParts
from generator.pieces import generatePiecesCutSurfaces
from generator.surfaces import generateSurfaces, generateSurfacesMeshes
from generator.tabs import generateTabsShells
from utils.utils import Configuration

import FreeCADTools.utils as FreeCADUtils

def generateElements(configuration: Configuration):
    def sliceAndGetSorteredResult(part, tool):
        sliceResult = FreeCADUtils.slicePart(part, [tool])
        
        if len(sliceResult) != 2:
            raise InvalidOperation
        
        if sliceResult[0].Volume > sliceResult[1].Volume:
            return sliceResult[0], sliceResult[1]
        else:
            return sliceResult[1], sliceResult[0]

    referenceSurface, topOffsetSurface, baseOffsetSurface = \
        generateSurfaces(configuration)
    
    document = FreeCADUtils.createNewDocument()

    print("Creando superficies")
    surfaceMesh, baseMesh = generateSurfacesMeshes(configuration, referenceSurface)

    body = FreeCADUtils.convertMeshToSolid(surfaceMesh)
    base = FreeCADUtils.convertMeshToSolid(baseMesh)

    # print("Creando tabs")
    # tabsShells = list(generateTabsShells(configuration, 
    #     referenceSurface, topOffsetSurface, baseOffsetSurface))

    # for tabShell, tabBaseShell in tabsShells:
    #     tabShellMesh = FreeCADUtils.convertMeshToSolid(tabShell)
    #     body, tab = sliceAndGetSorteredResult(body, tabShellMesh)

    #     tabShellMesh = FreeCADUtils.convertMeshToSolid(tabBaseShell)
    #     _, tabBase = sliceAndGetSorteredResult(base, tabShellMesh)
        
    #     tab = tab.fuse(tabBase)

    #     print("Tab creado")
    #     FreeCADUtils.addPartToDocument(tab)

    # print("Creando paredes y perforaciones")
    # walls, drills = generatePanelsParts(configuration, referenceSurface)
    # FreeCADUtils.addPartsToDocument(walls)

    # print("Perforando superficie")
    #for drill in drills:
    #    body = body.cut(drill)
    
    print("Creando planos y cortando superficie")
    cutSurfaces = list(generatePiecesCutSurfaces(configuration, referenceSurface))
    FreeCADUtils.addPartsToDocument(cutSurfaces)
    FreeCADUtils.addPartToDocument(body)
    #FreeCADUtils.addPartsToDocument(FreeCADUtils.slicePart(body, cutSurfaces))

    return document