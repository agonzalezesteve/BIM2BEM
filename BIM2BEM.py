
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

import ifcopenshell
import ifcopenshell.geom as geom
from skgeom import *
import skgeom
import pyclipper

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

SCALING_FACTOR = 1e3

def get_clipper(subjects, clips):
  clipper = pyclipper.Pyclipper()
  
  clipper.AddPaths(pyclipper.scale_to_clipper(subjects, SCALING_FACTOR), pyclipper.PT_SUBJECT, True)
  clipper.AddPaths(pyclipper.scale_to_clipper(clips, SCALING_FACTOR), pyclipper.PT_CLIP, True)
  
  return clipper

def execute_clipping(clipper, type):
  return pyclipper.scale_from_clipper(clipper.Execute(type, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD), SCALING_FACTOR)

def clipping(subjects, clips, operation):
  clipper = get_clipper(subjects, clips)
  
  if operation == "union": return execute_clipping(clipper, pyclipper.CT_UNION)
  if operation == "intersection": return execute_clipping(clipper, pyclipper.CT_INTERSECTION)
  if operation == "difference": return execute_clipping(clipper, pyclipper.CT_DIFFERENCE)
  
  union = execute_clipping(clipper, pyclipper.CT_UNION)
  intersection = execute_clipping(clipper, pyclipper.CT_INTERSECTION)
  
  if not intersection: return union
  
  clipper = get_clipper(union, intersection)
  
  return execute_clipping(clipper, pyclipper.CT_DIFFERENCE)

for plane, polygons in plane2polygons.items():
  solution = [polygons[0].coords.tolist()]
  for polygon in polygons[1:]:
    solution = clipping(solution, [polygon.coords.tolist()], "xor")
    
  print([plane.a(), plane.b(), plane.c(), plane.d()])
  for polygon in solution:
    print(list(map(lambda vertex: plane.to_3d(Point2(vertex[0], vertex[1])), polygon)))
  print("")
  # print([plane.a(), plane.b(), plane.c(), plane.d()])
  # for polygon in solution.polygons:
    # for vertex in polygon.outer_boundary().coords:
      # print(Point2(vertex[0], vertex[1]))
      # print(plane.to_3d(Point2(vertex[0], vertex[1])))
    # print("")
  # print("")
  
  plane2polygons[plane] = solution