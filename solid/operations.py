import numpy as np
from surfaces.iSurface import ISurface

from surfaces.meshSurface import MeshSurface


def getSurfaceVerticesGrid(surfaceF, xVect, yVect):
    vertices = []
    for x in xVect:
        row_vertices = []
        for y in yVect:
            row_vertices.append([x, y, surfaceF([x, y])])
        vertices.append(row_vertices)
    return vertices


def getPlaneTriangles(grid, ix, iy):
    return [[grid[ix][iy], grid[ix][iy+1], grid[ix+1][iy]],
            [grid[ix][iy+1], grid[ix+1][iy+1], grid[ix+1][iy]]]


def getSurfaceTriangleMesh(surfaceVerticesGrid):
    mesh_triangles = []
    for ix in range(len(surfaceVerticesGrid) - 1):
        for iy in range(len(surfaceVerticesGrid[ix]) - 1):
            mesh_triangles.extend(getPlaneTriangles(
                surfaceVerticesGrid, ix, iy))

    return mesh_triangles


def getLeadMesh(leadVertices1, leadVertices2, isCloseShape = False):

    def generateTriangle(i):
        return [[leadVertices1[i], leadVertices1[i+1], leadVertices2[i]] ,
                [leadVertices2[i], leadVertices1[i+1], leadVertices2[i+1]]]    

    if(len(leadVertices1) != len(leadVertices2)):
        raise IndexError

    mesh_triangles = []

    for i in range(-1 if isCloseShape else 0, len(leadVertices1)-1):
        mesh_triangles.extend(generateTriangle(i))

    return mesh_triangles


def getShellLeadsTriangleMesh(svg1, svg2):
    # svgX : surfaceVerticesGridX

    if(len(svg1) != len(svg2)):
        raise IndexError

    mesh_triangles = getLeadMesh(svg1[0],  svg2[0])

    mesh_triangles.extend(getLeadMesh(svg1[-1], svg2[-1]))

    lead_vertices1 = [v[0] for v in svg1]
    lead_vertices2 = [v[0] for v in svg2]

    mesh_triangles.extend(
        getLeadMesh(lead_vertices1, lead_vertices2))

    lead_vertices1 = [v[-1] for v in svg1]
    lead_vertices2 = [v[-1] for v in svg2]

    mesh_triangles.extend(
        getLeadMesh(lead_vertices1, lead_vertices2))

    return mesh_triangles

def getShellMeshFromSurfacesVertexs(sf1Vertexs, sf2Vertexs):
    mesh_triangles = []

    mesh_triangles.extend(getSurfaceTriangleMesh(sf1Vertexs))
    mesh_triangles.extend(getSurfaceTriangleMesh(sf2Vertexs))
    mesh_triangles.extend(getShellLeadsTriangleMesh(sf1Vertexs, sf2Vertexs))

    return mesh_triangles

def getShellMesh(surface1, surface2, xstart, xend, xsteps, ystart, yend, ysteps):

    xVect = np.linspace(xstart, xend, xsteps + 1)
    yVect = np.linspace(ystart, yend, ysteps + 1)    

    svg1 = getSurfaceVerticesGrid(surface1.F, xVect, yVect)
    svg2 = getSurfaceVerticesGrid(surface2.F, xVect, yVect)

    return getShellMeshFromSurfacesVertexs(svg1, svg2)

def getShellMesh(surface1 : MeshSurface, surface2 : MeshSurface):
    svg1, _ = surface1.getMeshData()
    svg2, _ = surface2.getMeshData()

    return getShellMeshFromSurfacesVertexs(svg1, svg2)

def surfaceToOffsetMeshSurface(surface: ISurface, offset, xVect, yVect):
    if surface is MeshSurface:
        raise TypeError
        
    surfaceF = lambda p: surface.FOffset(p, offset) if offset == 0 else surface.F

    surfaceVerticesGrid = getSurfaceVerticesGrid(surfaceF, xVect, yVect)
    surfaceTriangleMesh = getSurfaceTriangleMesh(surfaceVerticesGrid)

    return MeshSurface(surfaceVerticesGrid, surfaceTriangleMesh)

def surfaceToMeshSurface(surface: ISurface, xVect, yVect):
    return surfaceToOffsetMeshSurface(surface, 0, xVect, yVect)