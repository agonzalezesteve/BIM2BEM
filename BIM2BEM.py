
import ifcopenshell
import ifcopenshell.geom as geom
from skgeom import *
import skgeom
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

settings = ifcopenshell.geom.settings()
settings.set(settings.USE_WORLD_COORDS, True)

file = ifcopenshell.open('./test.ifc')

building_element2polygons = {}
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
    for key, value in plane2polygon_sets.items():
      is_opposite = key == plane.opposite()
      if key == plane or is_opposite:
        plane = key
        polygon_sets = value
        if is_opposite: points.reverse()
        break
    index = (0,1)[is_opposite]
    polygon_sets[index] = polygon_sets[index].join(skgeom.Polygon(list(map(lambda point: plane.to_2d(point), points))))

    plane2polygon_sets[plane] = polygon_sets
  
  building_element2polygons[building_element] = {plane: polygon_sets[0].symmetric_difference(polygon_sets[1]) for plane, polygon_sets in plane2polygon_sets.items()}

for building_element, plane2polygon_set in building_element2polygons.items():
  print(building_element)
  
  print(len(plane2polygon_set))
  for plane, polygon_set in plane2polygon_set.items():
    print(plane)
    for polygon in polygon_set.polygons: print(polygon)
    # for polygon in polygon_set.polygons: print(list(map(lambda vertex: plane.to_3d(Point2(vertex[0], vertex[1])), polygon.outer_boundary().coords)))

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