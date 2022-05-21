import numpy as np


def getSurfaceVerticesGrid(surface, xstart, xend, xsteps, ystart, yend, ysteps):
    vertices = []
    for x in np.linspace(xstart, xend, xsteps + 1):
        row_vertices = []
        for y in np.linspace(ystart, yend, ysteps + 1):
            row_vertices.append([x, y, surface.F([x, y])])
        vertices.append(row_vertices)
    return vertices


def getPlaneTriangles(grid, ix, iy):
    return [[grid[ix][iy], grid[ix][iy+1], grid[ix+1][iy]],
            [grid[ix][iy+1], grid[ix+1][iy+1], grid[ix+1][iy]]]


def gerSurfaceMesh(surfaceVerticesGrid):
    mesh_triangles = []
    for ix in range(len(surfaceVerticesGrid) - 1):
        for iy in range(len(surfaceVerticesGrid[ix]) - 1):
            mesh_triangles.extend(getPlaneTriangles(
                surfaceVerticesGrid, ix, iy))

    return mesh_triangles


def getLeadMesh(leadVertices1, leadVertices2):
    if(len(leadVertices1) != len(leadVertices2)):
        raise IndexError

    mesh_triangles = []

    for i in range(len(leadVertices1)-1):
        mesh_triangles.append(
            [leadVertices1[i], leadVertices1[i+1], leadVertices2[i]])
        mesh_triangles.append(
            [leadVertices2[i], leadVertices1[i+1], leadVertices2[i+1]])

    return mesh_triangles


def getShellLeadsMesh(svg1, svg2):
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


def getShellMesh(surface1, surface2, xstart, xend, xsteps, ystart, yend, ysteps):
    svg1 = getSurfaceVerticesGrid(
        surface1, xstart, xend, xsteps, ystart, yend, ysteps)
    svg2 = getSurfaceVerticesGrid(
        surface2, xstart, xend, xsteps, ystart, yend, ysteps)

    mesh_triangles = []

    mesh_triangles.extend(gerSurfaceMesh(svg1))
    mesh_triangles.extend(gerSurfaceMesh(svg2))
    mesh_triangles.extend(getShellLeadsMesh(svg1, svg2))

    return mesh_triangles
