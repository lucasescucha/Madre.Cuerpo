# Modificar el archivo freecad.pth a env/site-packages para que funcione la importación de FreeCAD
# pyright: reportMissingImports=false

import FreeCAD

import Mesh
import Part
import Draft

import BOPTools.SplitAPI as SplitAPI
from utils.utils import getRotationAngleAndAxis

TOLERANCE = 0.1

def convertToPolygon(polygon):
    vertices = [FreeCAD.Vector(v[0], v[1], v[2]) for v in polygon]
    wire = Part.makePolygon(vertices)
    return wire

def convertToFace(polygon):
    return Part.Face(polygon)

def toFreeCADVector(vector):
    return FreeCAD.Vector(vector[0], vector[1], vector[2])

def makeCompound(parts):
    return Part.makeCompound(parts)

def extrudePolygon(polygon, displacement):
    return polygon.extrude(toFreeCADVector(displacement))

def createPlane(length, width, pnt, dir):
    fcPnt = FreeCAD.Vector(pnt[0], pnt[1], pnt[2])
    fcDir = FreeCAD.Vector(dir[0], dir[1], dir[2])
    return Part.makePlane(length, width, fcPnt, fcDir)

def create3dText(text, size, extrusion):
    font = "/font/arial.ttf"
    text = Draft.makeShapeString(text, font, size)
    return extrudePolygon(text.Shape, extrusion)

def createText(text, pnt):
    fcPnt = FreeCAD.Vector(pnt[0], pnt[1], pnt[2])
    return Draft.makeText(text, fcPnt)

def slicePart(part, tools):
    result = SplitAPI.slice(part, tools, "Split")
    return result.Solids

def createCylinder(radius, lenght, center, direction):
    return Part.makeCylinder(radius, lenght, toFreeCADVector(center), 
            toFreeCADVector(direction), 360)

def createCone(radius1, radius2, lenght, center, direction):
    return Part.makeCone(radius1, radius2, lenght, toFreeCADVector(center), 
            toFreeCADVector(direction), 360)

def createSphere(radius, center):
    return Part.makeSphere(radius, toFreeCADVector(center))

def createBox(length, width, height, point, direction):
    return Part.makeBox(length,width,height, toFreeCADVector(point), 
            toFreeCADVector(direction))

def convertMeshToSolid(trianglesMesh):
    mesh = Mesh.Mesh(trianglesMesh)
    part = Part.Shape()
    part.makeShapeFromMesh(mesh.Topology, TOLERANCE)

    return Part.Solid(part)

def addPartsToDocument(parts):
    for part in parts:
        addPartToDocument(part)

def addMeshesToDocument(meshes):
    for mesh in meshes:
        Mesh.show(mesh)

def addMeshToDocument(mesh):
    Mesh.show(mesh)

def addPartToDocument(part):
    Part.show(part)


def createNewDocument():
    return FreeCAD.newDocument()

def saveDocument(document, path):
    document.saveAs(path)

def openDocument(path):
    return FreeCAD.open(path)