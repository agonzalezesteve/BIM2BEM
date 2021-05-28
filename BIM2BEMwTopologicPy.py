
import ifcopenshell
import ifcopenshell.geom as geom
from topologic import Vertex, Edge, Wire, Face, Shell, Cell, CellComplex, Cluster, Topology, Graph, Dictionary, Attribute, AttributeManager, VertexUtility, EdgeUtility, WireUtility, ShellUtility, CellUtility, TopologyUtility
import cppyy

def edgesByVertices(vertices):
  edges = []

  edges = cppyy.gbl.std.list[Edge.Ptr]()
  for i in range(len(vertices)-1):
    v1 = vertices[i]
    v2 = vertices[i+1]
    e1 = Edge.ByStartVertexEndVertex(v1, v2)
    edges.push_back(e1)

  # connect the last vertex to the first one
  v1 = vertices[len(vertices)-1]
  v2 = vertices[0]
  e1 = Edge.ByStartVertexEndVertex(v1, v2)
  edges.push_back(e1)

  return edges


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

settings = ifcopenshell.geom.settings()
settings.set(settings.USE_WORLD_COORDS, True)

ifc_file = ifcopenshell.open('./test.ifc')

build_elem_ccs = []
build_cc = None
for building_element in ifc_file.by_type('IfcBuildingElement'):
  if not (building_element.is_a('IfcWall') or building_element.is_a('IfcSlab')): continue

  shape = ifcopenshell.geom.create_shape(settings, building_element)
  geo = shape.geometry
  geo_vertices = geo.verts
  geo_faces = geo.faces

  topo_vertices = []
  for v in range(0, len(geo_vertices), 3):
    vertex = Vertex.ByCoordinates(geo_vertices[v], geo_vertices[v+1], geo_vertices[v+2])
    topo_vertices.append(vertex)

  topo_faces = cppyy.gbl.std.list[Face.Ptr]()
  for f in range(0, len(geo_faces), 3):
    face_vertices = []
    for v in geo_faces[f : f + 3]:
      vertex = topo_vertices[v]
      face_vertices.append(vertex)
    edges = edgesByVertices(face_vertices)
    face = Face.ByEdges(edges)
    topo_faces.push_back(face)
  
  cc = CellComplex.ByFaces(topo_faces, 0.0001)
  build_elem_ccs.append(cc)
  if build_cc is None:
    build_cc = cc
  else:
    build_cc = Topology.Merge(build_cc, cc)

ext_boundary = getSubTopologies(build_cc, CellComplex)[0].ExternalBoundary()
spaces = getSubTopologies(ext_boundary, Shell)
for i in range(len(spaces)):
  print(str(i+1)+". Space")
  space_faces = getSubTopologies(spaces[i], Face)
  for j in range(len(space_faces)):
    fVertices = getSubTopologies(space_faces[j], Vertex)
    if abs(fVertices[0].X() - fVertices[1].X()) < 1e-6 and abs(fVertices[0].X() - fVertices[2].X()) < 1e-6:
      print("  "+str(j+1)+". X: "+str(fVertices[0].X()))
    elif abs(fVertices[0].Y() - fVertices[1].Y()) < 1e-6 and abs(fVertices[0].Y() - fVertices[2].Y()) < 1e-6:
      print("  "+str(j+1)+". Y: "+str(fVertices[0].Y()))
    else:
      print("  "+str(j+1)+". Z: "+str(fVertices[0].Z()))