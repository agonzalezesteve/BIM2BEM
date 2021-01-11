
import ifcopenshell
import ifcopenshell.geom as geom
from skgeom import *
import skgeom
import pyclipper
import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

# poly_tree = pyclipper.PyPolyNode()
# poly_tree.Contour = [[40000, -0], [40000, 2500], [0, 2500]]
# poly_tree.IsOpen = True

# clipper = pyclipper.Pyclipper()

# clipper.AddPaths(pyclipper.PolyTreeToPaths(poly_tree), pyclipper.PT_SUBJECT, True)

# print(puta)

def are_equal(plane_a, plane_b):
  normal_a = plane_a.orthogonal_vector()
  normal_b = plane_b.orthogonal_vector()
  
  length_a = math.sqrt(normal_a.squared_length())
  length_b = math.sqrt(normal_b.squared_length())
  
  aux = normal_a * normal_b / length_a / length_b - 1
  if aux > 1e-6 or aux < -1e-6: return False
  
  aux = plane_a.d() / length_a - plane_b.d() / length_b
  if aux > 1e-6 or aux < -1e-6: return False
  
  return True

SCALING_FACTOR = 1e4

def get_poly_tree(plane, points):
  poly_tree = pyclipper.PyPolyNode()
  
  path = []
  for point in points:
    point = plane.to_2d(point)
    path.append([int(SCALING_FACTOR * point.x()), int(SCALING_FACTOR * point.y())])
  poly_tree.Contour = path
  poly_tree.IsOpen = True
  
  return poly_tree

def get_clipper(subject, clip):  
  clipper = pyclipper.Pyclipper()
  
  print(subject)
  
  clipper.AddPaths(subject, pyclipper.PT_SUBJECT, True)
  clipper.AddPaths(clip, pyclipper.PT_CLIP, True)
  
  return clipper

def execute_clipping(clipper, type):
  return clipper.Execute2(type, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD)

def clipping(subject, clip, boolean_operation):
  if not subject.Contour and not subject.Childs:
    if boolean_operation == "join": return clip
    if boolean_operation == "intersection": return subject
    if boolean_operation == "difference": return subject
    return clip
    
  if not clip.Contour and not clip.Childs:
    if boolean_operation == "join": return subject
    if boolean_operation == "intersection": return clip
    if boolean_operation == "difference": return subject
    return subject
  
  clipper = get_clipper(pyclipper.PolyTreeToPaths(subject), pyclipper.PolyTreeToPaths(clip))

  if boolean_operation == "join": return execute_clipping(clipper, pyclipper.CT_UNION)
  if boolean_operation == "intersection": return execute_clipping(clipper, pyclipper.CT_INTERSECTION)
  if boolean_operation == "difference": return execute_clipping(clipper, pyclipper.CT_DIFFERENCE)
  
  union = execute_clipping(clipper, pyclipper.CT_UNION)
  intersection = execute_clipping(clipper, pyclipper.CT_INTERSECTION)
  
  if not intersection: return union
  
  clipper = get_clipper(pyclipper.PolyTreeToPaths(union), pyclipper.PolyTreeToPaths(intersection))
  
  return execute_clipping(clipper, pyclipper.CT_DIFFERENCE)
  
settings = ifcopenshell.geom.settings()
settings.set(settings.USE_WORLD_COORDS, True)

file = ifcopenshell.open('./test.ifc')

building_element2poly_treee = {}
for building_element in file.by_type('IfcBuildingElement'):
  if not (building_element.is_a('IfcWall') or building_element.is_a('IfcSlab')): continue
  
  shape = ifcopenshell.geom.create_shape(settings, building_element)
  geo = shape.geometry
  vertices = geo.verts
  faces = geo.faces

  plane2poly_trees = {}
  for f in range(0, len(faces), 3):
    points = []
    for v in faces[f : f + 3]:
      points.append(skgeom.Point3(vertices[3*v], vertices[3*v+1], vertices[3*v+2]))
    plane = skgeom.Plane3(points[0], points[1], points[2])

    is_opposite = False
    poly_trees = [pyclipper.PyPolyNode(), pyclipper.PyPolyNode()]
    for building_element_plane, building_element_poly_trees in plane2poly_trees.items():
      is_opposite = are_equal(building_element_plane, plane.opposite())
      if is_opposite: points.reverse()
      if are_equal(building_element_plane, plane) or is_opposite:
        plane = building_element_plane
        poly_trees = building_element_poly_trees
        break
    index = (0,1)[is_opposite]
    poly_trees[index] = clipping(poly_trees[index], get_poly_tree(plane, points), "join")

    plane2poly_trees[plane] = poly_trees
  
  building_element2poly_treee[building_element] = {plane: clipping(poly_trees[0], poly_trees[1], "symmetric_difference") for plane, poly_trees in plane2poly_trees.items()}

for building_element, building_element_plane2poly_tree in building_element2poly_treee.items():
  print(building_element)
  for plane, poly_tree in building_element_plane2poly_tree.items():
    print(plane)
    pyclipper.PolyTreeToPaths(poly_tree)

print(puta)

def translate_vertex(vertex, old_plane, new_plane):
  return new_plane.to_2d(old_plane.to_3d(Point2(vertex[0], vertex[1])))
  
def translate_polygon_set(polygon_set, old_plane, new_plane, is_opposite):
  new_polygons = []
  
  for polygon in polygon_set.polygons:
    points = list(map(lambda vertex: translate_vertex(vertex, old_plane, new_plane), polygon.outer_boundary().coords))
    if is_opposite: points.reverse()
    outer = skgeom.Polygon(points)
    holes = []
    for hole in polygon.holes:
      points = list(map(lambda vertex: translate_vertex(vertex, old_plane, new_plane), hole.coords))
      if is_opposite: points.reverse()
      holes.append(skgeom.Polygon(points))
    new_polygons.append(skgeom.PolygonWithHoles(outer, holes))
  
  return skgeom.PolygonSet(new_polygons)
      
plane2polygon_set = {}
for building_element, building_element_plane2polygon_set in building_element2poly_treee.items():
  for plane, polygon_set in building_element_plane2polygon_set.items():
    for global_plane, global_polygon_set in plane2polygon_set.items():
      is_opposite = are_equal(global_plane, plane.opposite())
      if are_equal(global_plane, plane) or is_opposite:
        polygon_set = clipping(global_polygon_set, translate_polygon_set(polygon_set, plane, global_plane, is_opposite), "symmetric_difference")
        plane = global_plane
        break
        
    plane2polygon_set[plane] = polygon_set
    
count = 0
for plane, polygon_setpolygon_set in plane2polygon_set.items():
  print(count)
  print(plane)
  for polygon in polygon_set.polygons: 
    print("outer:")
    print(list(map(lambda vertex: plane.to_3d(Point2(vertex[0], vertex[1])), polygon.outer_boundary().coords)))
    # print(polygon.outer_boundary().coords)
    for hole in polygon.holes:
      print("hole:")
      print(list(map(lambda vertex: plane.to_3d(Point2(vertex[0], vertex[1])), hole.coords)))
      # print(hole.coords)
  print("")
  count += 1
  
print(puta)

def cross_product(u, v):
  return Vector3(Point3(u.z()*v.y(), u.x()*v.z(), u.y()*v.x()), Point3(u.y()*v.z(), u.z()*v.x(), u.x()*v.y()))

def do_intersect(polygon_i, polygon_j):
  for id_vertex_i, p2 in enumerate(polygon_i):
    p1 = polygon_i[id_vertex_i-1]
    p12 = p2 - p1
    for id_vertex_j, p4 in enumerate(polygon_j):
      p3 = polygon_j[id_vertex_j-1]
      p34 = p4 - p3

      if cross_product(p12, p34).squared_length() > 1e-6: continue

      p13 = p3 - p1

      if cross_product(p12, p13).squared_length() > 1e-6: continue

      p14 = p4 - p1

      k2 = p12 * p12
      k3 = p12 * p13
      k4 = p12 * p14

      intersection = [
        [0.0, p14.squared_length(), p12.squared_length()],
        [p13.squared_length(), p34.squared_length(), (p3 - p2).squared_length()],
        [p12.squared_length(), (p4 - p2).squared_length(), 0.0]
      ]

      if k3 < 1e-6:
        intersection_i = 0
      elif k3 < k2:
        intersection_i = 1
      else:
        intersection_i = 2

      if k4 < 1e-6:
        intersection_j = 0
      elif k4 < k2:
        intersection_j = 1
      else:
        intersection_j = 2

      if intersection[intersection_i][intersection_j] > 1e-6: return True
  return False

    # polygons.append(points)

  # row = []
  # col = []
  # data = []
  # for id_polygon_i, polygon_i in enumerate(polygons):
    # plane_i = skgeom.Plane3(polygon_i[0], polygon_i[1], polygon_i[2])
    # for id_polygon_j, polygon_j in enumerate(polygons[id_polygon_i+1:]):
      # plane_j = skgeom.Plane3(polygon_j[0], polygon_j[1], polygon_j[2])
      # if plane_i == plane_j and do_intersect(polygon_i, polygon_j):
        # row.append(id_polygon_i)
        # col.append(id_polygon_j+id_polygon_i+1)
        # data.append(1)
  # n_components, labels = connected_components(csgraph=csr_matrix((np.array(data), (np.array(row), np.array(col))), shape=(len(polygons), len(polygons))), directed=False, return_labels=True)
  # print(labels)

firsts = []
id, aux = 0, 0
for plane, polygons in building_element2polygons[element].items():
  solution = skgeom.PolygonSet()
  for polygon in polygons:
    if id == aux:
      print("a:")
      for a in solution.polygons:
        print(a)
      print("b:")
      print(polygon)
    solution = solution.symmetric_difference(polygon)
    if id == aux:
      print("c:")
      for c in solution.polygons:
        print(c)
      print(puta)
    id += 1

  for polygon in solution.polygons:
    firsts.append(list(map(lambda vertex: plane.to_3d(Point2(vertex[0], vertex[1])), polygon.outer_boundary().coords)))

print(puta)

for id_polygon_i, polygon_i in enumerate(firsts):
  print(id_polygon_i)
  for vertex in polygon_i:
    print(vertex)
  print("")

row = []
col = []
data = []
for id_polygon_i, polygon_i in enumerate(firsts):
  for id_polygon_j, polygon_j in enumerate(firsts[id_polygon_i+1:]):
    if do_intersect(polygon_i, polygon_j):
      print(str(id_polygon_i) + ": " + str(id_polygon_j+id_polygon_i+1))
      row.append(id_polygon_i)
      col.append(id_polygon_j+id_polygon_i+1)
      data.append(1)

n_components, labels = connected_components(csgraph=csr_matrix((np.array(data), (np.array(row), np.array(col))), shape=(len(firsts), len(firsts))), directed=False, return_labels=True)
# print(labels)