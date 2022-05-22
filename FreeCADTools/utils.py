# Modificar el archivo freecad.pth a env/site-packages para que funcione la importación de FreeCAD
# En .vscode settings.json se puede agregar "python. .analysis.extraPaths": ["<PATH TO FREECAD"] de
# forma de que pylance reconozca los paquetes

import FreeCAD

import Mesh
import Part
import Draft

import BOPTools.SplitAPI as SplitAPI

TOLERANCE = 0.1


def convertToPolygon(polygon, z):
    vertices = [FreeCAD.Vector(v[0], v[1], z) for v in polygon]
    wire = Part.makePolygon(vertices)
    return wire


def extrudePolygonInZ(polygon, zStart, zEnd):
    fcPolygon = convertToPolygon(polygon, zStart)
    return fcPolygon.extrude(FreeCAD.Vector(0, 0, zEnd - zStart))


def createPlane(length, width, pnt, dir):
    fcPnt = FreeCAD.Vector(pnt[0], pnt[1], pnt[2])
    fcDir = FreeCAD.Vector(dir[0], dir[1], dir[2])
    return Part.makePlane(length, width, fcPnt, fcDir)

def createText(text, pnt):
    fcPnt = FreeCAD.Vector(pnt[0], pnt[1], pnt[2])
    return Draft.makeText(text, fcPnt)

def slicePart(part, tools):
    result = SplitAPI.slice(part, tools, "Split")
    return result.Solids


def convertMeshToSolid(trianglesMesh):
    mesh = Mesh.Mesh(trianglesMesh)
    part = Part.Shape()
    part.makeShapeFromMesh(mesh.Topology, TOLERANCE)

    return Part.Solid(part)


def addPartsToDocument(parts):
    for part in parts:
        addPartToDocument(part)


def addPartToDocument(part):
    Part.show(part)


def createNewDocument():
    return FreeCAD.newDocument()


def saveDocument(document, path):
    document.saveAs(path)
