from fenics import *

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):
    # ylo boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool((near(x[0], 0, DOLFIN_EPS) or near(x[1], 0, DOLFIN_EPS)) and on_boundary)
    # Map yhi boundary (H) to ylo boundary (G)
    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1]
        y[2] = x[2]
        if near(x[0], 1, DOLFIN_EPS):
            y[0] = x[0] - 1.0
        if near(x[1], 1, DOLFIN_EPS):
            y[1] = x[1] - 1.0

# Create mesh and define function space
ngrids = 20
mesh = UnitCubeMesh.create(ngrids, ngrids, ngrids, CellType.Type_tetrahedron)
pds = PeriodicBoundary()
V = FunctionSpace(mesh, 'P', 1, constrained_domain=pds)

# Define boundary condition
class Back(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0, DOLFIN_EPS)

class Front(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1, DOLFIN_EPS)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0, DOLFIN_EPS)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 1, DOLFIN_EPS)

class BottomHeat(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] <= 0.5 and x[1] <= 0.5 and x[2] < DOLFIN_EPS and on_boundary)

class BottomNone(SubDomain):
    def inside(self, x, on_boundary):
        return bool((x[0] > 0.5 or x[1] > 0.5) and x[2] < DOLFIN_EPS and on_boundary)

class TopHeat(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] >= 0.5 and x[1] >= 0.5 and x[2] > (1-DOLFIN_EPS) and on_boundary)

class TopNone(SubDomain):
    def inside(self, x, on_boundary):
        return bool((x[0] < 0.5 or x[1] < 0.5) and x[2] > (1-DOLFIN_EPS) and on_boundary)

bk = Back()
ft = Front()
lt = Left()
rt = Right()
bh = BottomHeat()
bn = BottomNone()
th = TopHeat()
tn = TopNone()

bc0 = DirichletBC(V, Constant(0.0), bh)
bc1 = DirichletBC(V, Constant(1.0), th)
bcs = [bc0, bc1]

boundary_markers = MeshFunction("size_t", mesh, mesh.geometry().dim()-1)
boundary_markers.set_all(0)
bk.mark(boundary_markers, 1)
ft.mark(boundary_markers, 2)
lt.mark(boundary_markers, 3)
rt.mark(boundary_markers, 4)
bh.mark(boundary_markers, 5)
bn.mark(boundary_markers, 6)
th.mark(boundary_markers, 7)
tn.mark(boundary_markers, 8)

# Redefine boundary integration measure
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

# Define variational problem
kappa = 1.0
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0)
a = kappa*dot(grad(u), grad(v))*dx
L = f*v*dx

# Set linear solver parameters
linear_solver = 'Krylov'
prm = LinearVariationalSolver.default_parameters()
if linear_solver == 'Krylov':
    prm.linear_solver = 'gmres'
    prm.preconditioner = 'ilu'
    prm.krylov_solver.absolute_tolerance = 1e-10
    prm.krylov_solver.relative_tolerance = 1e-8
    prm.krylov_solver.maximum_iterations = 10000
else:
    prm.linear_solver = 'lu'

# Compute solution
u = Function(V)
solve(a == L, u, bcs, solver_parameters=prm)

# Evaluate integral of normal gradient over top boundary
n = FacetNormal(mesh)
total_flux = assemble(kappa*dot(grad(u), n)*ds)
print("total_flux = ", total_flux)
back       = assemble(kappa*dot(grad(u), n)*ds(1))
print("back       = ", back)
front      = assemble(kappa*dot(grad(u), n)*ds(2))
print("front      = ", front)
left       = assemble(kappa*dot(grad(u), n)*ds(3))
print("left       = ", left)
right      = assemble(kappa*dot(grad(u), n)*ds(4))
print("right      = ", right)
bottomheat = assemble(kappa*dot(grad(u), n)*ds(5))
print("bottomheat = ", bottomheat)
bottomnone = assemble(kappa*dot(grad(u), n)*ds(6))
print("bottomnone = ", bottomnone)
topheat    = assemble(kappa*dot(grad(u), n)*ds(7))
print("topheat    = ", topheat)
topnone    = assemble(kappa*dot(grad(u), n)*ds(8))
print("topnone    = ", topnone)

# Save solution to file in VTK format
vtkfile = File('periodic_3d/solution.pvd')
vtkfile << u
