"""
Test problem is chosen to give an exact solution at all nodes of the mesh.
  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary
  u_D = 0, 0<x<0.5 && y==0;
  u_D = 1, 0.5<x<1 && y==1;
    f = 0
"""

from fenics import *
import matplotlib.pyplot as plt

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)
    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 1.0
        y[1] = x[1]

# Create mesh and define function space
mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'P', 1, constrained_domain=PeriodicBoundary())

# Define boundary condition
class BottomLeft(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] <= 0.5 and x[1] < DOLFIN_EPS and on_boundary)

class BottomRight(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] > 0.5 and x[1] < DOLFIN_EPS and on_boundary)

class TopLeft(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] <= 0.5 and x[1] > (1-DOLFIN_EPS) and on_boundary)

class TopRight(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] > 0.5 and x[1] > (1-DOLFIN_EPS) and on_boundary)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0, DOLFIN_EPS)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1, DOLFIN_EPS)

bl = BottomLeft()
br = BottomRight()
tl = TopLeft()
tr = TopRight()
ll = Left()
rr = Right()

bc0 = DirichletBC(V, Constant(0.0), bl)
bc1 = DirichletBC(V, Constant(1.0), tr)
bcs = [bc0, bc1]

boundary_markers = MeshFunction("size_t", mesh, mesh.geometry().dim()-1)
boundary_markers.set_all(0)
bl.mark(boundary_markers, 1)
br.mark(boundary_markers, 2)
tl.mark(boundary_markers, 3)
tr.mark(boundary_markers, 4)
ll.mark(boundary_markers, 5)
rr.mark(boundary_markers, 6)

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
total_flux   = assemble(kappa*dot(grad(u), n)*ds)
print("total_flux   = ", total_flux)
bottom_left  = assemble(kappa*dot(grad(u), n)*ds(1))
print("bottom_left  = ", bottom_left)
bottom_right = assemble(kappa*dot(grad(u), n)*ds(2))
print("bottom_right = ", bottom_right)
top_left     = assemble(kappa*dot(grad(u), n)*ds(3))
print("top_left     = ", top_left)
top_right    = assemble(kappa*dot(grad(u), n)*ds(4))
print("top_right    = ", top_right)
left         = assemble(kappa*dot(grad(u), n)*ds(5))
print("left         = ", left)
right        = assemble(kappa*dot(grad(u), n)*ds(6))
print("right        = ", right)

# Plot solution and mesh
plot(u)
plot(mesh)

# Save solution to file in VTK format
vtkfile = File('periodic/solution.pvd')
vtkfile << u

# Hold plot
plt.show()