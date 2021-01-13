
import ifcopenshell
import ifcopenshell.geom as geom
from skgeom import *
import skgeom
import pyclipper
import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def are_equal(plane_a, plane_b):
  normal_a = plane_a.orthogonal_vector()
  normal_b = plane_b.orthogonal_vector()
  
  length_a = math.sqrt(normal_a.squared_length())
  length_b = math.sqrt(normal_b.squared_length())
  
  if not abs(float(normal_a * normal_b / length_a / length_b - 1)) < 1e-6: return False
  
  if not abs(float(plane_a.d() / length_a - plane_b.d() / length_b)) < 1e-6: return False

  return True

SCALING_FACTOR = 1e4

def get_poly_tree(plane, points):
  poly_tree = pyclipper.PyPolyNode()
  
  polygon2 = skgeom.Polygon(list(map(lambda point: plane.to_2d(point), points)))
  poly_tree.Contour = list(map(lambda coord: [int(SCALING_FACTOR * coord[0]), int(SCALING_FACTOR * coord[1])] , polygon2.coords))
  poly_tree.IsOpen = True
  
  return poly_tree

def get_clipper(subject, clip):  
  clipper = pyclipper.Pyclipper()
    
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
  
  if not intersection.Contour and not intersection.Childs:
    return union
  
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

# for building_element, building_element_plane2poly_tree in building_element2poly_treee.items():
  # print(building_element)
  # for plane, poly_tree in building_element_plane2poly_tree.items():
    # print(plane)
    # for polygon2 in pyclipper.PolyTreeToPaths(poly_tree):
      # print(list(map(lambda coord: plane.to_3d(skgeom.Point2(coord[0] / SCALING_FACTOR, coord[1] / SCALING_FACTOR)), polygon2)))

# print(puta)

def translate_poly_tree(poly_tree, parent, old_plane, new_plane, is_opposite):
  points = list(map(lambda coord: old_plane.to_3d(skgeom.Point2(coord[0] / SCALING_FACTOR, coord[1] / SCALING_FACTOR)), poly_tree.Contour))
  if is_opposite: points.reverse()
  new_poly_tree = get_poly_tree(new_plane, points)
  if parent: new_poly_tree.Parent = parent
  
  new_poly_tree.IsHole = poly_tree.IsHole
  new_poly_tree.IsOpen = poly_tree.IsOpen
  
  depths = [0]
  for child in poly_tree.Childs:
    new_child = translate_poly_tree(child, new_poly_tree, old_plane, new_plane, is_opposite)
    new_poly_tree.Childs.append(new_child)
    depths.append(new_child.depth + 1)
  new_poly_tree.depth = max(depths)
  
  return new_poly_tree
      
plane2poly_tree = {}
for building_element, building_element_plane2poly_tree in building_element2poly_treee.items():
  for plane, poly_tree in building_element_plane2poly_tree.items():
    for global_plane, global_poly_tree in plane2poly_tree.items():
      is_opposite = are_equal(global_plane, plane.opposite())
      if are_equal(global_plane, plane) or is_opposite:
        aux = translate_poly_tree(poly_tree, None, plane, global_plane, is_opposite)
        poly_tree = clipping(global_poly_tree, translate_poly_tree(poly_tree, None, plane, global_plane, is_opposite), "symmetric_difference")
        plane = global_plane
        break
        
    plane2poly_tree[plane] = poly_tree
    
# count = 0
# for plane, poly_tree in plane2poly_tree.items():
  # print(count)
  # print(plane)
  # for polygon2 in pyclipper.PolyTreeToPaths(poly_tree):
    # print(list(map(lambda coord: plane.to_3d(skgeom.Point2(coord[0] / SCALING_FACTOR, coord[1] / SCALING_FACTOR)), polygon2)))
  # print("")
  # count += 1

def path2polygon3(path, plane):
  return list(map(lambda coord: plane.to_3d(skgeom.Point2(coord[0] / SCALING_FACTOR, coord[1] / SCALING_FACTOR)), path))

def add_polygon3s(poly_tree, plane, firsts):
  contour = poly_tree.Contour
  if not contour:
    for child in poly_tree.Childs: add_polygon3s(child, plane, firsts)
  else:
    outer_boundary = path2polygon3(contour, plane)
    holes = list(map(lambda hole: path2polygon3(hole.Contour, plane), poly_tree.Childs))
    
    firsts.append([plane, outer_boundary, holes])
    
    for hole in poly_tree.Childs: 
      for child in hole.Childs: add_polygon3s(child, plane, firsts)

firsts = []
for plane, poly_tree in plane2poly_tree.items():
  add_polygon3s(poly_tree, plane, firsts)
  
for polygon3 in firsts:
  print(polygon3)
  
def has_on_plane(plane, point):
  normal = plane.orthogonal_vector()
  distance = (normal * (point - skgeom.Point3(0, 0, 0)) + plane.d()) / math.sqrt(normal.squared_length())
  
  return abs(float(distance)) < 1e-3

def get_segment3s(paths, plane):
  # print(paths)
  segment3s = []
  
  for path in paths:
    prev_point3 = path[-1]
    has_prev_point3 = has_on_plane(plane, prev_point3)
    for point3 in path:
      if has_on_plane(plane, point3):
        if has_prev_point3: segment3s.append(skgeom.Segment3(point3, prev_point3))
        has_prev_point3 = True
      else:
        has_prev_point3 = False
      prev_point3 = point3
  
  return segment3s

def segment2_is_vertical(segment2):
  u = segment2.to_vector()
  v = skgeom.Vector2(0, 1)
  
  return abs(float(u.x()*v.y()-u.y()*v.x())) < 1e-3
  
def polygons3_do_intersect(polygon3_i, polygon3_j):      
  segment3s_i = get_segment3s(([polygon3_i[1]] + polygon3_i[2]), polygon3_j[0])
  if not segment3s_i: return False
  plane_i = polygon3_i[0]
  segment3s_j = get_segment3s(([polygon3_j[1]] + polygon3_j[2]), plane_i)
  if not segment3s_j: return False
  
  for segment3_i in segment3s_i:
    segment2_i = skgeom.Segment2(plane_i.to_2d(segment3_i.source()), plane_i.to_2d(segment3_i.target()))
    is_vertical = segment2_is_vertical(segment2_i)
    
    point2_min_i, point2_max_i = segment2_i.min(), segment2_i.max()
    min_i = point2_min_i.x() if not is_vertical else point2_min_i.y()
    max_i = point2_max_i.x() if not is_vertical else point2_max_i.y()
    for segment3_j in segment3s_j:
      segment2_j = skgeom.Segment2(plane_i.to_2d(segment3_j.source()), plane_i.to_2d(segment3_j.target()))
      point2_min_j, point2_max_j = segment2_j.min(), segment2_j.max()
      min_j = point2_min_j.x() if not is_vertical else point2_min_j.y()
      max_j = point2_max_j.x() if not is_vertical else point2_max_j.y()

      if max_j > min_i and max_i > min_j: return True
  
  return False


row = []
col = []
data = []
for id_polygon_i, polygon3_i in enumerate(firsts):
  for id_polygon_j, polygon3_j in enumerate(firsts[id_polygon_i+1:]):
    if not polygons3_do_intersect(polygon3_i, polygon3_j): continue
    
    row.append(id_polygon_i)
    col.append(id_polygon_j+id_polygon_i+1)
    data.append(1)
    
n_components, labels = connected_components(csgraph=csr_matrix((np.array(data), (np.array(row), np.array(col))), shape=(len(firsts), len(firsts))), directed=False, return_labels=True)

print(labels)

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