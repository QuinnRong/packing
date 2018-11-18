from fenics import *

num = 1

mesh = UnitSquareMesh.create(num, num, CellType.Type_triangle)			# 三角形
File("output/triangle.pvd") << mesh

mesh = UnitSquareMesh.create(num, num, CellType.Type_quadrilateral)		# 四边形
File("output/quadrilateral.pvd") << mesh

mesh = UnitCubeMesh.create(num, num, num, CellType.Type_tetrahedron)	# 四面体
File("output/tetrahedron.pvd") << mesh

mesh = UnitCubeMesh.create(num, num, num, CellType.Type_hexahedron)		# 六面体
File("output/hexahedron.pvd") << mesh
