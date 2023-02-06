import gmsh
import sys

gmsh.initialize()
gmsh.model.add("task_0_kubik")


mesh_param = 0.01
coordinates_center = [0,1,0]
radius = 1
point = []
curves = []
array_curves=[]
A = [[-1,-1,-1],[-1,1,-1],[1,1,-1],[1,-1,-1],[-1,-1,1],[-1,1,1],[1,1,1],[1,-1,1]]

for i in range(8):
    point.append(gmsh.model.geo.addPoint(coordinates_center[0] + A[i][0]*radius, coordinates_center[1] + A[i][1]*radius, coordinates_center[2] + A[i][2]*radius, mesh_param))

for i in range(4):
    gmsh.model.geo.addLine(point[i], point[i+4])
    gmsh.model.geo.addLine(point[i%4], point[(i+1)%4])
    gmsh.model.geo.addLine(point[i%4+4], point[(i+1)%4+4])

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(3)

gmsh.write("cube.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()