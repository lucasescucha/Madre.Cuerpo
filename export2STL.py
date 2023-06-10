import os
import shutil

import FreeCADTools.utils as FreeCADUtils

document = FreeCADUtils.openDocument("output/output.FCStd")

if os.path.exists("output/stl"):
    shutil.rmtree("output/stl", ignore_errors=True)

os.makedirs("output/stl")

for obj in document.Objects:
    print("Exportando objeto", obj.Label)
    filename = "output/stl/" + obj.Label + ".stl"
    obj.Shape.exportStl(filename)    

input("Presione una tecla para terminar...")