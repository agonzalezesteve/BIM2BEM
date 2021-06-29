
import bpy
import bmesh
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
import ifcopenshell
import ifcopenshell.geom as geom
import topologic
from topologic import Vertex, Edge, Wire, Face, Shell, Cell, CellComplex, Cluster, Topology, Graph, Dictionary, Attribute, AttributeManager, VertexUtility, EdgeUtility, WireUtility, FaceUtility, ShellUtility, CellUtility, TopologyUtility
import cppyy
from blenderbim.bim.ifc import IfcStore

def classByType(argument):
  switcher = {
    1: Vertex,
    2: Edge,
    4: Wire,
    8: Face,
    16: Shell,
    32: Cell,
    64: CellComplex,
    128: Cluster }

  return switcher.get(argument, Topology)

def fixTopologyClass(topology):
  topology.__class__ = classByType(topology.GetType())

  return topology

def getSubTopologies(topology, subTopologyClass):
  pointer = subTopologyClass.Ptr
  values = cppyy.gbl.std.list[pointer]()
  if subTopologyClass == Vertex:
    _ = topology.Vertices(values)
  elif subTopologyClass == Edge:
    _ = topology.Edges(values)
  elif subTopologyClass == Wire:
    _ = topology.Wires(values)
  elif subTopologyClass == Face:
    _ = topology.Faces(values)
  elif subTopologyClass == Shell:
    _ = topology.Shells(values)
  elif subTopologyClass == Cell:
    _ = topology.Cells(values)
  elif subTopologyClass == CellComplex:
    _ = topology.CellComplexes(values)

  py_list = []
  i  =  values.begin()
  while (i != values.end()):
    py_list.append(fixTopologyClass(Topology.DeepCopy(i.__deref__())))
    _ = i.__preinc__()

  return py_list

def vertexIndex(v, vertices, tolerance):
  index = None

  v._class__ = Vertex
  for i in range(len(vertices)):
    vertices[i].__class__ = Vertex
    d = VertexUtility.Distance(v, vertices[i])
    if d <= tolerance:
      index = i
      break

  return index

def meshData(topology):
  type = classByType(topology.GetType())
  if type == Cluster or type == CellComplex or type == Cell or type == Shell:
    topFaces = getSubTopologies(topology, Face)
    topVertices = getSubTopologies(topology, Vertex)
  else:
    topFaces = []
    topVertices = []

  vertices = []
  for aVertex in topVertices:
    vertices.append(tuple([aVertex.X(), aVertex.Y(), aVertex.Z()]))

  faces = []
  for aFace in topFaces:
    wires = getSubTopologies(aFace, Wire)
    wire = wires[0]
    faceVertices = getSubTopologies(wire, Vertex)
    tempList = []
    for aVertex in faceVertices:
      index = vertexIndex(aVertex, topVertices, 1e-4)
      tempList.append(index)
    faces.append(tuple(tempList))

  return [vertices, faces]

settings = ifcopenshell.geom.settings()
settings.set(settings.USE_BREP_DATA, True)
settings.set(settings.SEW_SHELLS, True)
settings.set(settings.USE_WORLD_COORDS, True)

ifc_file = IfcStore.get_file()

build_elem_cells = []
build_top = None
for build_elem in ifc_file.by_type('IfcBuildingElement'):
  if not (build_elem.is_a('IfcWall') or build_elem.is_a('IfcSlab')): continue

  shape = ifcopenshell.geom.create_shape(settings, build_elem)
  brepString = shape.geometry.brep_data
  build_elem_top = fixTopologyClass(Topology.ByString(brepString))
  faces = cppyy.gbl.std.list[Face.Ptr]()
  for face in getSubTopologies(build_elem_top, Face):
    faces.push_back(face)
  build_elem_cell = Cell.ByFaces(faces, 1e-4)
  build_elem_cells.append(build_elem_cell)

  int_bounds = cppyy.gbl.std.list[Shell.Ptr]()
  build_elem_cell.InternalBoundaries(int_bounds)
  for shell in int_bounds:
    md = meshData(shell)
    mesh = bpy.data.meshes.new(name="IfcConnectionSurfaceGeometry")
    mesh.from_pydata(md[0], [], md[1])
    obj = object_data_add(bpy.context, mesh)
    
  if build_top is None:
    build_top = build_elem_cell
  else:
    build_top = Topology.Merge(build_top, build_elem_cell)

int_bounds = cppyy.gbl.std.list[Face.Ptr]()
getSubTopologies(build_top, CellComplex)[0].InternalBoundaries(int_bounds)
for face in int_bounds:
  faces = cppyy.gbl.std.list[Face.Ptr]()
  faces.push_back(face)
  shell = Shell.ByFaces(faces, 1e-4)
  md = meshData(shell)
  mesh = bpy.data.meshes.new(name="IfcConnectionSurfaceGeometry")
  mesh.from_pydata(md[0], [], md[1])
  obj = object_data_add(bpy.context, mesh)
  cells = cppyy.gbl.std.list[Cell.Ptr]()
  FaceUtility.AdjacentCells(face, build_top, cells)

# ext_boundary = getSubTopologies(build_cc, CellComplex)[0].ExternalBoundary()
# spaces = getSubTopologies(ext_boundary, Shell)
# print(len(spaces))
# int_boundaries = cppyy.gbl.std.list[Face.Ptr]()
# getSubTopologies(build_cc, CellComplex)[0].InternalBoundaries(int_boundaries)
# for j, face in enumerate(list(int_boundaries)):
  # fVertices = getSubTopologies(face, Vertex)
  # if abs(fVertices[0].X() - fVertices[1].X()) < 1e-6 and abs(fVertices[0].X() - fVertices[2].X()) < 1e-6:
    # print("  "+str(j+1)+". X: "+str(fVertices[0].X()))
  # elif abs(fVertices[0].Y() - fVertices[1].Y()) < 1e-6 and abs(fVertices[0].Y() - fVertices[2].Y()) < 1e-6:
    # print("  "+str(j+1)+". Y: "+str(fVertices[0].Y()))
  # else:
    # print("  "+str(j+1)+". Z: "+str(fVertices[0].Z()))