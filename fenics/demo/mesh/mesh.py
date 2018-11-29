from fenics import *

ngrids = 2

types = {1:{"dim":2,"num_per_cell":2,"type":CellType.Type_triangle},
         2:{"dim":2,"num_per_cell":1,"type":CellType.Type_quadrilateral},
         3:{"dim":3,"num_per_cell":6,"type":CellType.Type_tetrahedron},
         4:{"dim":3,"num_per_cell":1,"type":CellType.Type_hexahedron}}

for idx in types:
    if types[idx]["dim"] == 2:
        mesh = UnitSquareMesh.create(ngrids, ngrids, types[idx]["type"])
    if types[idx]["dim"] == 3:
        mesh = UnitCubeMesh.create(ngrids, ngrids, ngrids, types[idx]["type"])
    File("output/mesh_"+str(idx)+".pvd") << mesh

    materials = MeshFunction("double", mesh, mesh.geometry().dim())
    for (i, cell) in enumerate(cells(mesh)):
        materials[i] = i / types[idx]["num_per_cell"]
    File("output/material_"+str(idx)+".pvd") << materials
