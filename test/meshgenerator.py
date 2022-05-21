from solver.offsetSurfaceSolver import OffsetSurfaceSolver
from surfaces.motherSurface import MotherSurface
from surfaces.planeSurface import PlaneSurface
from surfaces.offset.motherOffsetSurface import MotherOffsetSurface
from surfaces.shiftDecorator import ShiftDecorator

import matplotlib.pyplot as plt
import numpy as np

surface = ShiftDecorator(MotherSurface(0.01, 0.01), [0, 10, 0])

offset = MotherOffsetSurface(surface, 5)

solver = OffsetSurfaceSolver(surface)

x_axis = range(-10, 100)

surface = [surface.F(np.array([x, 0])) for x in x_axis]
offset1 = [offset.F(np.array([x, 0])) for x in x_axis]
offset2 = [solver.getOffsetSurfaceZ(np.array([x, 0]), 10)
           for x in x_axis]
offset3 = [solver.getOffsetSurfaceZ(np.array([x, 0]), 15)
           for x in x_axis]

fig = plt.figure()
ax = fig.add_subplot()

ax.plot(x_axis, surface, color='tab:blue')
ax.plot(x_axis, offset1, color='tab:orange')
ax.plot(x_axis, offset2, color='tab:green')
ax.plot(x_axis, offset3, color='tab:red')

plt.xlim(-10, 100)
plt.ylim(-10, 100)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
