import os

from generator.generator import generateElements

import FreeCADTools.utils as FreeCADUtils

from configuration import configurationData
from utils.utils import Configuration


def run(checkParts=False):

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
    
    input("Presione una tecla para terminar...")

run(False)