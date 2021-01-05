
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

import ifcopenshell
import ifcopenshell.geom as geom
from skgeom import *
import skgeom
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def draw(p):
	fig = plt.figure()

	faces = []
	for facet in p.facets:
		face = []
		first = facet.halfedge()
		point = first.vertex().point()
		face.append([float(point.x()), float(point.y()), float(point.z())])
		he = first.next()
		while he is not first:
			point = he.vertex().point()
			face.append([float(point.x()), float(point.y()), float(point.z())])
			he = he.next()
		faces.append(face)

	arr = np.asarray(faces)
	ax = Axes3D(fig)
	ax.add_collection3d(mplot3d.art3d.Poly3DCollection(faces, alpha=0.3, edgecolor='red'))

	all_points = arr.reshape(-1, 9)
	xlim = np.min(all_points[:, 0]), np.max(all_points[:, 0])
	ylim = np.min(all_points[:, 1]), np.max(all_points[:, 1])
	zlim = np.min(all_points[:, 2]), np.max(all_points[:, 2])
	lims = np.asarray([xlim, ylim, zlim])
	lims = np.min(lims[:, 0]), np.max(lims[:, 1])
	ax.set_xlim3d(*lims)
	ax.set_ylim3d(*lims)
	ax.set_zlim3d(*lims)
	plt.show()

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

settings = ifcopenshell.geom.settings()
settings.set(settings.USE_WORLD_COORDS, True)

file = ifcopenshell.open('./test.ifc')
  
plane2polygons = {}
building_element2polyhedron_3 = {}
for building_element in file.by_type('IfcBuildingElement'):
  if not (building_element.is_a('IfcWall') or building_element.is_a('IfcSlab')): continue
  
  shape = ifcopenshell.geom.create_shape(settings, building_element)
  geo = shape.geometry
  vertices = geo.verts
  faces = geo.faces

  polygons = []
  for f in range(0, len(faces), 3):
    points = []
    for v in faces[f : f + 3]:
      points.append(skgeom.Point3(vertices[3*v], vertices[3*v+1], vertices[3*v+2]))
    polygons.append(points)
    
  row = []
  col = []
  data = []
  for id_polygon_i, polygon_i in enumerate(polygons):
    plane_i = skgeom.Plane3(polygon_i[0], polygon_i[1], polygon_i[2])
    for id_polygon_j, polygon_j in enumerate(polygons[id_polygon_i+1:]):
      plane_j = skgeom.Plane3(polygon_j[0], polygon_j[1], polygon_j[2])
      if plane_i == plane_j and do_intersect(polygon_i, polygon_j):
        row.append(id_polygon_i)
        col.append(id_polygon_j+id_polygon_i+1)
        data.append(1)
  n_components, labels = connected_components(csgraph=csr_matrix((np.array(data), (np.array(row), np.array(col))), shape=(len(polygons), len(polygons))), directed=False, return_labels=True)
  print(labels)
        
        
    # plane = skgeom.Plane3(points[0], points[1], points[2])

    # polygons = []
    # for key, value in plane2polygons.items():
      # is_opposite = key == plane.opposite
      # if key == plane or is_opposite:
        # plane = key
        # polygons = value
        # if is_opposite: points.reverse()
        # break
    # polygons.append(skgeom.Polygon(list(map(lambda point: plane.to_2d(point), points))))

    # plane2polygons[plane] = polygons
  
  # building_element2polyhedron_3[building_element] = polyhedron

print(puta)

firsts = []
id, aux = 0, 3
for plane, polygons in plane2polygons.items():
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