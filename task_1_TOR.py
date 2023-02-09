import gmsh
import sys

gmsh.initialize()

tor_radiuses = [10,3,6]
coord_center = [0,0,0]
dim = 3

gmsh.model.add("task_1_9K330")

TORS = [gmsh.model.occ.addTorus(coord_center[0], coord_center[1], coord_center[2], tor_radiuses[0], tor_radiuses[i]) for i in [1,2]]
gmsh.model.occ.cut([(dim, TORS[1])], [(dim, TORS[0])], dim)
gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.3)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("task_1_9K330.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
