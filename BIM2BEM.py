
import math
import ifcopenshell
import ifcopenshell.geom as geom
from skgeom import *
import skgeom
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

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


def append_points(obj, points):
  if isinstance(obj, skgeom._skgeom.Polygon):
    for coord in obj.coords: points.append(skgeom.Point2(coord[0], coord[1]))
      
    return True
    
  if isinstance(obj, skgeom._skgeom.PolygonWithHoles):
    append_points(obj.outer_boundary(), points)
    for hole in obj.holes: append_points(hole, points)
      
    return True
    
  if isinstance(obj, skgeom._skgeom.PolygonSet):
    for polygon in obj.polygons: append_points(polygon, points)
    
    return True

def replace_points(obj, grid):
  if isinstance(obj, skgeom._skgeom.Polygon):
    points = []
    
    for coord in obj.coords:
      point = skgeom.Point2(coord[0], coord[1])
      points.append([x for x in grid if skgeom.squared_distance(point, x) < 1e-6][0])
    
    return skgeom.Polygon(points)
    
  if isinstance(obj, skgeom._skgeom.PolygonWithHoles):
    outer_boundary = replace_points(obj.outer_boundary(), grid)
    holes = list(map(lambda hole: replace_points(hole, grid), obj.holes))
    
    return skgeom.PolygonWithHoles(outer_boundary, holes)
    
  if isinstance(obj, skgeom._skgeom.PolygonSet):
    polygons = list(map(lambda polygon: replace_points(polygon, grid), obj.polygons))
    
    return skgeom.PolygonSet(polygons)

def remove_collinear_vertices(obj):
  if isinstance(obj, skgeom._skgeom.Polygon):
    points = []
    
    for coord in obj.coords:
      point = skgeom.Point2(coord[0], coord[1])
      if len(points) > 0 and skgeom.squared_distance(points[-1], point) < 1e-6: points.pop()
      elif len(points) > 1 and skgeom.collinear(points[-2], points[-1], point): points.pop()
      points.append(point)
    
    if len(points) > 0 and skgeom.squared_distance(points[-1], points[0]) < 1e-6: points.pop()
    elif len(points) > 1 and skgeom.collinear(points[-2], points[-1], points[0]): points.pop()
    
    return skgeom.Polygon(points)
    
  if isinstance(obj, skgeom._skgeom.PolygonWithHoles):
    outer_boundary = remove_collinear_vertices(obj.outer_boundary())
    holes = list(map(lambda hole: remove_collinear_vertices(hole), obj.holes))
    holes = [hole for hole in holes if not hole.area() > -1e-6]
    
    return skgeom.PolygonWithHoles(outer_boundary, holes)
    
  if isinstance(obj, skgeom._skgeom.PolygonSet):
    polygons = list(map(lambda polygon: remove_collinear_vertices(polygon), obj.polygons))

    return skgeom.PolygonSet(polygons)

def get_area(obj):
  if isinstance(obj, skgeom._skgeom.Polygon):
    return obj.area()
    
  if isinstance(obj, skgeom._skgeom.PolygonWithHoles):
    area = get_area(obj.outer_boundary())
    for hole in obj.holes: area += get_area(hole)
    
    return area
    
  if isinstance(obj, skgeom._skgeom.PolygonSet):
    area = 0.0
    for polygon in obj.polygons: area += get_area(polygon)

    return area
    
def print_polygon(obj):
  if isinstance(obj, skgeom._skgeom.Polygon):
    print(obj.coords)
    
  if isinstance(obj, skgeom._skgeom.PolygonWithHoles):
    print("outer:")
    print_polygon(obj.outer_boundary())
    print("holes:")
    for hole in obj.holes: print_polygon(hole)
    
  if isinstance(obj, skgeom._skgeom.PolygonSet):
    for polygon in obj.polygons: print_polygon(polygon)

    return area
    
def clipping(subjects, clips, boolean_operation):  
  points = []
  
  append_points(subjects, points)
  append_points(clips, points)
  
  row = []
  col = []
  data = []
  for id_i, point_i in enumerate(points):
    for id_j, point_j in enumerate(points[id_i+1:]):
      if skgeom.squared_distance(point_i, point_j) < 1e-6:
        row.append(id_i)
        col.append(id_j+id_i+1)
        data.append(1)
  n_components, labels = connected_components(csgraph=csr_matrix((np.array(data), (np.array(row), np.array(col))), shape=(len(points), len(points))), directed=False, return_labels=True)
  
  new_points = []
  for x in range(0, n_components):
    indices = (np.where(labels == x))[0]
    if len(indices) == 1:
      new_points.append(points[indices[0]])
    else:
      old_points = list(map(lambda id: points[id], indices))
      new_points.append(skgeom.centroid(skgeom.Polygon(old_points)))
  
  subjects = replace_points(subjects, points)
  clips = replace_points(clips, points)

  method_to_call = getattr(subjects, boolean_operation)
  result = method_to_call(clips)

  # return result
  return remove_collinear_vertices(result)
  
settings = ifcopenshell.geom.settings()
settings.set(settings.USE_WORLD_COORDS, True)

file = ifcopenshell.open('./test.ifc')

building_element2polygon_set = {}
for building_element in file.by_type('IfcBuildingElement'):
  if not (building_element.is_a('IfcWall') or building_element.is_a('IfcSlab')): continue
  
  shape = ifcopenshell.geom.create_shape(settings, building_element)
  geo = shape.geometry
  vertices = geo.verts
  faces = geo.faces

  plane2polygon_sets = {}
  for f in range(0, len(faces), 3):
    points = []
    for v in faces[f : f + 3]:
      points.append(skgeom.Point3(vertices[3*v], vertices[3*v+1], vertices[3*v+2]))
    plane = skgeom.Plane3(points[0], points[1], points[2])

    is_opposite = False
    polygon_sets = [skgeom.PolygonSet(), skgeom.PolygonSet()]
    for building_element_plane, building_element_polygon_sets in plane2polygon_sets.items():
      is_opposite = are_equal(building_element_plane, plane.opposite())
      if is_opposite: points.reverse()
      if are_equal(building_element_plane, plane) or is_opposite:
        plane = building_element_plane
        polygon_sets = building_element_polygon_sets
        break
    index = (0,1)[is_opposite]
    polygon_sets[index] = clipping(polygon_sets[index], skgeom.Polygon(list(map(lambda point: plane.to_2d(point), points))), "join")

    plane2polygon_sets[plane] = polygon_sets
  
  building_element2polygon_set[building_element] = {plane: clipping(polygon_sets[0], polygon_sets[1], "symmetric_difference") for plane, polygon_sets in plane2polygon_sets.items()}

# for building_element, building_element_plane2polygon_set in building_element2polygon_set.items():
  # print(building_element)
  # for plane, polygon_set in building_element_plane2polygon_set.items():
    # print(plane)
    # for polygon in polygon_set.polygons: 
      # print("outer:")
      # print(list(map(lambda vertex: plane.to_3d(Point2(vertex[0], vertex[1])), polygon.outer_boundary().coords)))
      # for hole in polygon.holes:
        # print("hole:")
        # print(list(map(lambda vertex: plane.to_3d(Point2(vertex[0], vertex[1])), hole.coords)))

# print(puta)

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
for building_element, building_element_plane2polygon_set in building_element2polygon_set.items():
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