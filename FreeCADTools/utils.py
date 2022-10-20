# Modificar el archivo freecad.pth a env/site-packages para que funcione la importación de FreeCAD
# pyright: reportMissingImports=false

import FreeCAD

import Mesh
import Part
import Draft

import BOPTools.SplitAPI as SplitAPI

TOLERANCE = 0.1

def convertToPolygon(polygon):
    vertices = [FreeCAD.Vector(v[0], v[1], v[2]) for v in polygon]
    wire = Part.makePolygon(vertices)
    return wire

def toFreeCADVector(vector):
    return FreeCAD.Vector(vector[0], vector[1], vector[2])

def extrudePolygon(polygon, displacement):
    return polygon.extrude(toFreeCADVector(displacement))

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

def createCylinder(radius, lenght, center, direction):
    return FreeCAD.makeCylinder(radius, lenght, toFreeCADVector(center), toFreeCADVector(direction), 360)

def convertMeshToSolid(trianglesMesh):
    mesh = Mesh.Mesh(trianglesMesh)
    part = Part.Shape()
    part.makeShapeFromMesh(mesh.Topology, TOLERANCE)

    return Part.Solid(part)

def createMesh(trianglesMesh):
    return Mesh.Mesh(trianglesMesh)

def addPartsToDocument(parts):
    for part in parts:
        addPartToDocument(part)

def addMeshesToDocument(meshes):
    for mesh in meshes:
        Mesh.show(mesh)


def addPartToDocument(part):
    Part.show(part)


def createNewDocument():
    return FreeCAD.newDocument()


def saveDocument(document, path):
    document.saveAs(path)
