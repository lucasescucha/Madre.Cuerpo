from msilib.schema import Error
from sqlite3 import InternalError
from generator.surfaces import generateSurfaces
from generator.utils import calculatePiecesData
from utils.utils import Configuration
from configuration import configurationData

def checkConfiguration():
    configuration = Configuration(configurationData)

    referenceSurface, _, _ = generateSurfaces(configuration)

    surfDims = configuration.surface.dimensions
    moldConfig = configuration.manufacture.mold
    puzzleConfig = configuration.manufacture.mold.puzzle

    _, _, xPieces,_ = calculatePiecesData(referenceSurface, 
        surfDims, puzzleConfig.pieces, "x")

    _, _, yPieces,_ = calculatePiecesData(referenceSurface, 
        surfDims, puzzleConfig.pieces, "y")

    print("Piezas en x:", xPieces, "/Piezas en y:", yPieces)

    pnWidth = moldConfig.panels.dimensions.width
    pnHeight = moldConfig.panels.dimensions.height

    print("Piezas por panel en x", pnWidth, "Piezas por panel en y", pnHeight)

    if xPieces%pnWidth != 0 or yPieces%pnHeight != 0:
        raise InternalError

checkConfiguration()