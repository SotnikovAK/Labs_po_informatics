import vtk
import numpy as np
import gmsh
import math
import os


class CalcMesh:

    def __init__(self, nodes_coords, tetrs_points, velocities):
        self.nodes = np.array([nodes_coords[0::3], nodes_coords[1::3], nodes_coords[2::3]])
        self.norm_raspred = np.linalg.norm(self.nodes, axis=0)
        self.velocity = (self.norm_raspred * (velocities[1] - velocities[0]) / np.max(self.norm_raspred) + velocities[0]) * self.nodes / self.norm_raspred

        self.a = - (0.5 * np.linalg.norm(self.velocity, axis=0) ** 1.5)

        self.tetrs = np.array([tetrs_points[0::4], tetrs_points[1::4], tetrs_points[2::4], tetrs_points[3::4]])
        self.tetrs -= 1

    def move(self, tau):
        self.velocity += self.a * tau
        self.nodes += self.velocity * tau

    def snapshot(self, snap_number):
        unstructuredGrid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()

        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(3)
        vel.SetName("velocity")

        for i in range(0, len(self.nodes[0])):
            points.InsertNextPoint(self.nodes[0, i], self.nodes[1, i], self.nodes[2, i])
            vel.InsertNextTuple((self.velocity[0, i], self.velocity[1, i], self.velocity[2, i]))

        unstructuredGrid.SetPoints(points)
        unstructuredGrid.GetPointData().AddArray(vel)

        for i in range(0, len(self.tetrs[0])):
            tetr = vtk.vtkTetra()
            for j in range(0, 4):
                tetr.GetPointIds().SetId(j, self.tetrs[j, i])
            unstructuredGrid.InsertNextCell(tetr.GetCellType(), tetr.GetPointIds())

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(unstructuredGrid)
        writer.SetFileName("CAT-" + str(snap_number) + ".vtu")
        writer.Write()


gmsh.initialize()
try:
    path = os.path.dirname(os.path.abspath(__file__))
    gmsh.merge(os.path.join(path, 'geometric cat.stl'))
except:
    print("Could not load STL mesh: bye!")
    gmsh.finalize()
    exit(-1)

angle = 40
forceParametrizablePatches = False
includeBoundary = True
curveAngle = 180
gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary, forceParametrizablePatches,
                                 curveAngle * math.pi / 180.)
gmsh.model.mesh.createGeometry()

s = gmsh.model.getEntities(2)
l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
gmsh.model.geo.addVolume([l])

gmsh.model.geo.synchronize()

f = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(f, "F", "4")
gmsh.model.mesh.field.setAsBackgroundMesh(f)

gmsh.model.mesh.generate(3)

nodeTags, nodesCoord, parametricCoord = gmsh.model.mesh.getNodes()

GMSH_TETR_CODE = 4
tetrsNodesTags = None
elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements()
for i in range(0, len(elementTypes)):
    if elementTypes[i] != GMSH_TETR_CODE:
        continue
    tetrsNodesTags = elementNodeTags[i]

if tetrsNodesTags is None:
    print("Can not find tetra data. Exiting.")
    gmsh.finalize()
    exit(-2)

print("The model has %d nodes and %d tetrs" % (len(nodeTags), len(tetrsNodesTags) / 4))

for i in range(0, len(nodeTags)):
    assert (i == nodeTags[i] - 1)
assert (len(tetrsNodesTags) % 4 == 0)

temperatures = np.array([1, 500])
velocities = np.array([6, 2])
mesh = CalcMesh(nodesCoord, tetrsNodesTags, velocities)
for i in range(10):
    mesh.snapshot(i)
    mesh.move(0.05)

gmsh.finalize()