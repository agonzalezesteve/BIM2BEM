
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

import ifcopenshell
import ifcopenshell.geom as geom
from skgeom import *
import skgeom

points = []
points.append(skgeom.Point3(0, 0, 1))
p3 = skgeom.Point3(0, 1, 1)
points.append(p3)
points.append(skgeom.Point3(1, 0, 1))
plane = skgeom.Plane3(points[0], points[1], points[2])

print(p3)
print(plane.has_on(p3))
p2 = plane.to_2d(p3)
print(p2)
print(plane.to_3d(p2))

print(puta)

settings = ifcopenshell.geom.settings()
settings.set(settings.USE_WORLD_COORDS, True)

file = ifcopenshell.open('./test.ifc')

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

plane2polygons = {}
building_element2polyhedron_3 = {}
for building_element in file.by_type('IfcBuildingElement'):
  if not (building_element.is_a('IfcWall') or building_element.is_a('IfcSlab')): continue
  
  polyhedron = skgeom.Polyhedron_3()
  
  shape = ifcopenshell.geom.create_shape(settings, building_element)
  geo = shape.geometry
  
  vertices = geo.verts
  faces = geo.faces
  for f in range(0, len(faces), 3):
    points = []
    for v in faces[f : f + 3]:
      points.append(skgeom.Point3(vertices[3*v], vertices[3*v+1], vertices[3*v+2]))
    polyhedron.make_triangle(points[0], points[1], points[2])
    plane = skgeom.Plane3(points[0], points[1], points[2])
    
    polygons = []
    for key, value in plane2polygons.items():
      if key == plane or key == plane.opposite: 
        plane = key
        polygons = value
        break
    polygons.append(skgeom.Polygon(list(map(lambda point: plane.to_2d(point), points))))
    
    plane2polygons[plane] = polygons
  
  building_element2polyhedron_3[building_element] = polyhedron

for plane, polygons in plane2polygons.items():
  solution = skgeom.PolygonSet([polygons[0]])
  for polygon in polygons[1:]:
    solution = solution.symmetric_difference(polygon)
  
  print([plane.a(), plane.b(), plane.c(), plane.d()])
  for polygon in solution.polygons:
    for vertex in polygon.outer_boundary().vertices:
      print(plane.to_3d(vertex))
  print("")
  
  plane2polygons[plane] = solution

print(puta)

settings = ifcopenshell.geom.settings()
settings.set(settings.USE_WORLD_COORDS, True)

file = ifcopenshell.open('./test.ifc')

def three2two(inv_tf, point3d):
  vertex = inv_tf * point3d
  
  return [vertex.x(), vertex.y()]

def two2three(tf, point2d):
  vertex = tf * openstudio.Point3d(point2d[0], point2d[1], 0.0)
  
  return [vertex.x(), vertex.y(), vertex.z()]

SCALING_FACTOR = 1e4

def get_clipper(subjects2D, clips2D):
  clipper = pyclipper.Pyclipper()
  
  clipper.AddPaths(pyclipper.scale_to_clipper(subjects2D, SCALING_FACTOR), pyclipper.PT_SUBJECT, True)
  clipper.AddPaths(pyclipper.scale_to_clipper(clips2D, SCALING_FACTOR), pyclipper.PT_CLIP, True)
  
  return clipper

def union2D(clipper):
  return pyclipper.scale_from_clipper(clipper.Execute(pyclipper.CT_UNION, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD), SCALING_FACTOR)

def intersection2D(clipper):
  return pyclipper.scale_from_clipper(clipper.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD), SCALING_FACTOR)

def difference2D(clipper):
  return pyclipper.scale_from_clipper(clipper.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD), SCALING_FACTOR)

def clipping2D(subjects2D, clips2D, operation):
  clipper = get_clipper(subjects2D, clips2D)
  
  if operation == "union": return union2D(clipper)
  if operation == "intersection": return intersection2D(clipper)
  if operation == "difference": return difference2D(clipper)
  
  union = union2D(clipper)
  intersection = intersection2D(clipper)
  
  if not intersection: return union
  
  clipper = get_clipper(union, intersection)
  
  return difference2D(clipper)

def clipping3D(subjects3D, clips3D, operation):
  tf = openstudio.Transformation.alignFace(clips3D[0])
  inv_tf = tf.inverse()
  
  subjects2D = list(map(lambda polygon: list(map(lambda vertex: three2two(inv_tf, vertex), polygon)), subjects3D))
  clips2D = list(map(lambda polygon: list(map(lambda vertex: three2two(inv_tf, vertex), polygon)), clips3D))
  
  solution = clipping2D(subjects2D, clips2D, operation)

  return list(map(lambda polygon: list(map(lambda vertex: two2three(tf, point2d), polygon)), solution))

plane2polygons = {}
building_elements = file.by_type('IfcBuildingElement')
for building_element in building_elements:
  if not (building_element.is_a('IfcWall') or building_element.is_a('IfcSlab')): continue
  
  shape = ifcopenshell.geom.create_shape(settings, building_element)
  geo = shape.geometry
    
  vertices = geo.verts
  faces = geo.faces
  for f in range(0, len(faces), 3):
    point3ds = []
    for v in faces[f : f + 3]:
      point3ds.append(openstudio.Point3d(vertices[3*v], vertices[3*v+1], vertices[3*v+2]))
    plane = openstudio.Plane(point3ds)
    
    polygons = []
    for key, value in plane2polygons.items():
      if key.equal(plane) or key.reverseEqual(plane): 
        plane = key
        polygons = value
        break
    polygons.append(point3ds)
    
    plane2polygons[plane] = polygons

plane2inv_tfs = {plane: openstudio.Transformation.alignFace(polygons[0]).inverse() for plane, polygons in plane2polygons.items()}

for plane, polygons in plane2polygons.items():
  plane2polygons[plane] = list(map(lambda polygon: list(map(lambda vertex: three2two(plane2inv_tfs[plane], vertex), polygon)), polygons))
  
for plane, polygons in plane2polygons.items():
  solution = [polygons[0]]
  for polygon in polygons[1:]:
    solution = clipping2D(solution, [polygon], "xor")
  
  print([plane.a(), plane.b(), plane.c(), plane.d()])
  for polygon in solution:
    print(list(map(lambda vertex: two2three(plane2inv_tfs[plane].inverse(), vertex), polygon)))
  print("")
  
  plane2polygons[plane] = solution