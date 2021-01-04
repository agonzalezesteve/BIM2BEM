
import ifcopenshell
import ifcopenshell.geom as geom
import pyclipper
import openstudio

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