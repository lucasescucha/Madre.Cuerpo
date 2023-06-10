import os
import shutil

import FreeCADTools.utils as FreeCADUtils
import Mesh

document = FreeCADUtils.openDocument("output/output.FCStd")

if os.path.exists("output/stl"):
    shutil.rmtree("output/stl")

os.makedirs("output/stl")

for obj in document.Objects:
    print("Exportando objeto", obj.Label)
    filename = "output/stl/" + obj.Label + ".stl"
    Mesh.export([obj], filename)

input("Presione una tecla para terminar...")