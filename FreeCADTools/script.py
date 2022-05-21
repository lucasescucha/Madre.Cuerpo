import Part
import Mesh
import FreeCAD

import json

import os
CURR_DIR = os.path.dirname(os.path.realpath(__file__))
print(CURR_DIR)
print(os.listdir("../"))

mesh = []

with open("C:/Users/Lucas/OneDrive/Documentos/Arte/Madre/CAD/Conjunto encastrado/MeshGenerator/mesh.json") as f_in:
    mesh = json.load(f_in)

part = Part.Shape()
part.makeShapeFromMesh(Mesh.Mesh(mesh).Topology, 1)

solid = Part.Solid(part)

Part.show(solid)