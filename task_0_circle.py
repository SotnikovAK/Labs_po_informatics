import gmsh
import sys

gmsh.initialize()

gmsh.model.add("task_0_circle")
mesh_param = 0.01
coordinates_center = [0,0,0]
radius = 1
point = []
curves = []
array_curves=[]
A = [[0,0],[1,0],[0,1],[-1,0],[0,-1]]
B = [[0,0],[1,2],[2,3],[3,4],[4,1]]
for i in range(5):
    point.append(gmsh.model.geo.addPoint(coordinates_center[0] + A[i][0]*radius, coordinates_center[1]+ A[i][1]*radius, coordinates_center[2], mesh_param))

for i in range(1,5):
    curves.append(gmsh.model.geo.addCircleArc(point[B[i][0]], point[0], point[B[i][1]]))
for i in range(4):
    array_curves.append(curves[i])

gmsh.model.geo.addCurveLoop(array_curves, 1)


gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("task_0_circle.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()