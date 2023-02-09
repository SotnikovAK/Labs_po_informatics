import gmsh
import math
import os
import sys

gmsh.initialize()
gmsh.write('eurofighter_typhoon.msh')

gmsh.merge(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'eurofighter_typhoon.stl'))

a = 800
b = 500
ugol = 180
gmsh.model.mesh.classifySurfaces(ugol * math.pi/180., False, False , 0)
gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop([gmsh.model.getEntities(2)[i][1] for i in range(len(gmsh.model.getEntities(2)))])])

A = [[a,0],[a,b],[0,b],[0,0]]
point = []
for i in range(len(A)):
    point.append(gmsh.model.geo.addPoint( A[i][0], A[i][1], -50, 0.01))

for i in range(len(A)):
    gmsh.model.geo.addLine(point[i], point[(i+1)%len(A)])
    

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
