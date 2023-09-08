import FreeCADTools.utils as FreeCADUtils
import Part

import numpy as np

document = FreeCADUtils.createNewDocument()

r = 1
lenght = 5
center = np.array([2,3,-4])

direction = np.array([1,1,0])
direction = direction/np.linalg.norm(direction)

t = np.array([-1,1,0])
t = t/np.linalg.norm(t)

n = np.cross(direction, t)
n = n/np.linalg.norm(n)

def createHousing(position, radius, lenght, direction, normal, tangent, closed):
    c1 = FreeCADUtils.createCylinder(radius, lenght, position, direction)
    b1 = FreeCADUtils.createBox(2*radius, radius, lenght, position + tangent*radius + normal*radius, direction)

    body = c1.fuse(b1)

    if closed:
        s1 = FreeCADUtils.createSphere(radius, position + direction*lenght)
        c2 = FreeCADUtils.createCylinder(radius, radius, position + direction*lenght, normal)

        body = body.fuse(s1).fuse(c2)

    return body
    

Part.show(createHousing(center, r, 0.5*r, lenght, direction, n, t, False))
Part.show(createHousing(center, r, 0.5*r, lenght, direction, n, t, True))

FreeCADUtils.saveDocument(document, "output/output.FCStd")
