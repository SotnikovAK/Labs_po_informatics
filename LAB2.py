import vtk
import numpy as np
import gmsh
import math
import os


class CalcMesh:
    def __init__(self, nodes_coords, tetrs_points):
        self.nodes = np.array([nodes_coords[0::3], nodes_coords[1::3], nodes_coords[2::3]])
        self.smth = np.power(self.nodes[0, :], 2) + np.power(self.nodes[1, :], 2)

        self.velocity = np.zeros(shape=(3, int(len(nodes_coords) / 3)), dtype=np.double)
        self.x_0, self.y_0 = 0, 0
        for i in range(0, len(self.nodes[0])):
            self.x_0 += self.nodes[0][i]
            self.y_0 += self.nodes[1][i]
        self.x_0 /= len(self.nodes[0])
        self.y_0 /= len(self.nodes[0])

        self.modul = np.power(np.power(self.nodes[0] - self.x_0, 2) + np.power(self.nodes[1] - self.y_0, 2), 1/2)
        self.alpha = []
        for i in range(0, len(self.nodes[0])):
            if self.nodes[0][i] - self.x_0 > 0:
                self.alpha.append(np.arctan((self.nodes[1][i] - self.y_0) / (self.nodes[0][i] - self.x_0)))
            else:
                self.alpha.append(np.pi + np.arctan((self.nodes[1][i] - self.y_0) / (self.nodes[0][i] - self.x_0)))


        self.tetrs = np.array([tetrs_points[0::4], tetrs_points[1::4], tetrs_points[2::4], tetrs_points[3::4]])
        self.tetrs -= 1

        self.velocity[0] = -self.modul * np.sin(self.alpha)
        self.velocity[1] = self.modul * np.cos(self.alpha)
        self.velocity[2] = 0

    def move(self):
        self.nodes[0] = self.modul * np.cos(self.alpha) + self.x_0
        self.nodes[1] = self.modul * np.sin(self.alpha) + self.y_0
        self.velocity[0] = - self.modul * np.sin(self.alpha)
        self.velocity[1] = self.modul * np.cos(self.alpha)
        for i in range(0, len(self.alpha)):
            self.alpha[i] += 0.1

        self.smth = np.power(self.nodes[0, :], 2) + np.power(self.nodes[1, :], 2)

    def snapshot(self, snap_number):
        unstructuredGrid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()

        smth = vtk.vtkDoubleArray()
        smth.SetName("smth")

        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(3)
        vel.SetName("vel")


        for i in range(0, len(self.nodes[0])):
            points.InsertNextPoint(self.nodes[0, i], self.nodes[1, i], self.nodes[2, i])
            smth.InsertNextValue(self.smth[i])
            vel.InsertNextTuple((self.velocity[0, i], self.velocity[1, i], self.velocity[2, i]))

        unstructuredGrid.SetPoints(points)

        unstructuredGrid.GetPointData().AddArray(smth)
        unstructuredGrid.GetPointData().AddArray(vel)

        for i in range(0, len(self.tetrs[0])):
            tetr = vtk.vtkTetra()
            for j in range(0, 4):
                tetr.GetPointIds().SetId(j, self.tetrs[j, i])
            unstructuredGrid.InsertNextCell(tetr.GetCellType(), tetr.GetPointIds())

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(unstructuredGrid)
        writer.SetFileName("cat_step" + str(snap_number) + ".vtu")
        writer.Write()


gmsh.initialize()

try:
    path = os.path.dirname(os.path.abspath(__file__))
    gmsh.merge(os.path.join(path, 'geometric cat.stl'))
except:
    print("Could not load STL mesh: bye!")
    gmsh.finalize()
    exit(-1)

angle =10
forceParametrizablePatches = True
includeBoundary = False
curveAngle = 180
gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary, forceParametrizablePatches,
                                 curveAngle * math.pi / 180.)
gmsh.model.mesh.createGeometry()

s = gmsh.model.getEntities(2)
l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
gmsh.model.geo.addVolume([l])

gmsh.model.geo.synchronize()

f = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(f, "F", "1")
gmsh.model.mesh.createGeometry()
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

mesh = CalcMesh(nodesCoord, tetrsNodesTags)
mesh.snapshot(0)

gmsh.finalize()
for i in range(1, 100):
    mesh.move()
    mesh.snapshot(i)
