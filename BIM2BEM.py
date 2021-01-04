
import ifcopenshell
import ifcopenshell.geom as geom
import pyclipper
import openstudio

settings = ifcopenshell.geom.settings()
settings.set(settings.USE_WORLD_COORDS, True)

file = ifcopenshell.open('path/test.ifc')

def three2two(inv_trans, point3d):
  vertex = inv_trans * point3d
  return [vertex.x(), vertex.y()]

def two2three(trans, point2d):
  vertex = trans * openstudio.Point3d(point2d[0], point2d[1], 0.0)
  return [vertex.x(), vertex.y(), vertex.z()]

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
      are_reversed = key.reverseEqual(plane)
      if key.equal(plane) or are_reversed: 
        plane = key
        polygons = value
        if are_reversed: point3ds.reverse()
        break
    polygons.append(point3ds)
    
    plane2polygons[plane] = polygons

plane2inv_trans = {plane: openstudio.Transformation.alignFace(polygons[0]).inverse() for plane, polygons in plane2polygons.items()}

for plane, polygons in plane2polygons.items():
  plane2polygons[plane] = list(map(lambda polygon: list(map(lambda vertex: three2two(plane2inv_trans[plane], vertex), polygon)), polygons))
  
for plane, polygons in plane2polygons.items():  
  solution = [polygons[0]]
  for polygon in polygons[1:]:
    clipper = pyclipper.Pyclipper()
    clipper.AddPaths(pyclipper.scale_to_clipper(solution), pyclipper.PT_SUBJECT, True)
    clipper.AddPath(pyclipper.scale_to_clipper(polygon), pyclipper.PT_CLIP, True)
    solution = pyclipper.scale_from_clipper(clipper.Execute(pyclipper.CT_XOR, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD))
  
  print([plane.a(), plane.b(), plane.c(), plane.d()])
  for polygon in solution:
    print(list(map(lambda vertex: two2three(plane2inv_trans[plane].inverse(), vertex), polygon)))
  print("")
  
  plane2polygons[plane] = solution