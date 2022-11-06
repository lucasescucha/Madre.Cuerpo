import numpy as np

from solid.operations import getShellMesh, surfaceToMeshSurface, surfaceToOffsetMeshSurface
from surfaces.iSurface import ISurface
from surfaces.motherSurface import MotherSurface
from surfaces.offset.motherOffsetSurface import MotherOffsetSurface
from surfaces.shiftDecorator import ShiftDecorator
from surfaces.utils import calculateSamplingPoints
from utils.utils import Configuration


def generateSurfaces(configuration: Configuration):
    def calculateReferenceSurface():
        surfaceConfig = configuration.surface
        zOffset = surfaceConfig.shape.zeroHeight

        # z = zOffset + a*x^2 -  b*y^2

        # y=0 -> z = zOffset + a*x^2 -> height = zOffset + a*(width/2)^2
        # a = (height - zOffset) / (width/2)^2
        a = (surfaceConfig.dimensions.height - zOffset) / \
            ((surfaceConfig.dimensions.width/2)**2)

        # x=0 -> z = zOffset - b*y^2 -> 0 = zOffset - b*(depth/2)^2
        # b = zOffset / (depth/2)^2
        b = zOffset / ((surfaceConfig.dimensions.depth/2)**2)

        return ShiftDecorator(MotherSurface(a, b),  np.array([0, 0, zOffset]))

    thickness = configuration.manufacture.mold.puzzle.thickness
    tabsBaseThickness = configuration.manufacture.mold.puzzle.tabs.baseThickness

    baseOffset = thickness + tabsBaseThickness

    referenceSurface = calculateReferenceSurface()
    topOffsetSurface = MotherOffsetSurface(referenceSurface, thickness)
    baseOffsetSurface = MotherOffsetSurface(referenceSurface, baseOffset)

    return referenceSurface, topOffsetSurface, baseOffsetSurface

def generateSurfacesMeshes(configuration: Configuration, referenceSurface: ISurface):

    surface = configuration.surface.dimensions
    grid = configuration.manufacture.grid

    thickness = configuration.manufacture.mold.puzzle.thickness
    tabsBaseThickness = configuration.manufacture.mold.puzzle.tabs.baseThickness

    baseOffset = thickness + tabsBaseThickness

    xstart, xend = 0, surface.width/2
    xVect = calculateSamplingPoints(xstart, xend, referenceSurface, "x", grid.width)

    ystart, yend = (-surface.depth/2), (surface.depth / 2)
    yVect = calculateSamplingPoints(ystart, yend, referenceSurface, "y", grid.height)

    meshBottomSurface = surfaceToMeshSurface(referenceSurface, xVect, yVect)
    meshTopSurface = surfaceToOffsetMeshSurface(referenceSurface, thickness, xVect, yVect)
    meshBaseSurface = surfaceToOffsetMeshSurface(referenceSurface, baseOffset, xVect, yVect)

    surfaceMesh = getShellMesh(meshBottomSurface, meshTopSurface)
    baseMesh = getShellMesh(meshTopSurface, meshBaseSurface)

    return surfaceMesh, baseMesh

