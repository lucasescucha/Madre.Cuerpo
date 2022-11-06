import os
import shutil
from generator.generator import generateElements

import FreeCADTools.utils as FreeCADUtils

from configuration import configurationData
from utils.utils import Configuration


def run(checkParts=False, savePiecesStl=False):

    configuration = Configuration(configurationData)

    if not configuration.check():
        raise ArithmeticError

    document = generateElements(configuration)

    if not os.path.exists("output"):
        os.makedirs("output")
        
    if os.path.exists("output/output.FCStd"):
        os.remove("output/output.FCStd")

    if checkParts:
        if not document.Objects.all(lambda p: p.isValid()):
            raise SystemError

    FreeCADUtils.saveDocument(document, "output/output.FCStd")

    if savePiecesStl:
        if os.path.exists("output/stl"):
            shutil.rmtree("output/stl", ignore_errors=True)

        os.makedirs("output/stl")

        for obj in document.Objects:
            filename = "output/stl/" + obj.Label + ".stl"
            obj.Shape.exportStl(filename)

run(False, False)