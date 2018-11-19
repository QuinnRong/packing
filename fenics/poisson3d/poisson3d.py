from fenics import *
import numpy as np
import time

# global parameters
dim = 100
TC = [1, 10]

nscale = 2
ngrids = dim*nscale

linear_solver = 'Krylov'
err_rt = 1e-6
err_at = 1e-8

types = {1:{"dim":2,"num_per_cell":2,"type":CellType.Type_triangle},
         2:{"dim":2,"num_per_cell":1,"type":CellType.Type_quadrilateral},
         3:{"dim":3,"num_per_cell":6,"type":CellType.Type_tetrahedron},
         4:{"dim":3,"num_per_cell":1,"type":CellType.Type_hexahedron}}
cell_type = types[4]

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):
    # xlo and ylo boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool((near(x[0], 0, DOLFIN_EPS) or near(x[1], 0, DOLFIN_EPS)) and on_boundary)
    # Map xhi and yhi boundary (H) to ylo boundary (G)
    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1]
        y[2] = x[2]
        if near(x[0], 1, DOLFIN_EPS):
            y[0] = x[0] - 1.0
        if near(x[1], 1, DOLFIN_EPS):
            y[1] = x[1] - 1.0

# Create mesh and define function space
mesh = UnitCubeMesh.create(ngrids, ngrids, ngrids, cell_type["type"])
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

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 0, DOLFIN_EPS)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 1, DOLFIN_EPS)

bk = Back()
ft = Front()
lt = Left()
rt = Right()
bm = Bottom()
tp = Top()

bc0 = DirichletBC(V, Constant(0.0), bm)
bc1 = DirichletBC(V, Constant(1.0), tp)
bcs = [bc0, bc1]

boundary_markers = MeshFunction("size_t", mesh, mesh.geometry().dim()-1)
boundary_markers.set_all(0)
bk.mark(boundary_markers, 1)
ft.mark(boundary_markers, 2)
lt.mark(boundary_markers, 3)
rt.mark(boundary_markers, 4)
bm.mark(boundary_markers, 5)
tp.mark(boundary_markers, 6)

# Redefine boundary integration measure
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

def GetPrm():
    # Set linear solver parameters
    prm = LinearVariationalSolver.default_parameters()
    if linear_solver == 'Krylov':
        prm.linear_solver = 'gmres'
        prm.preconditioner = 'ilu'
        prm.krylov_solver.absolute_tolerance = err_at
        prm.krylov_solver.relative_tolerance = err_rt
        prm.krylov_solver.maximum_iterations = 100000
    else:
        prm.linear_solver = 'lu'
    return prm

def GetKappa(idx):
    dirname = "./" + str(int(idx))
    kappa_input = []
    for j in range(0, dim): 
        filename = dirname+"/3D_"+str(int(idx))+"_"+str(int(j))+".dat"
        kappa_layer = np.loadtxt(filename, unpack=True).ravel()
        kappa_input = np.append(kappa_input, kappa_layer)

    volume_fraction = 0
    for index in range(kappa_input.size):
        if (kappa_input[index] < 0.5):
            kappa_input[index] = TC[0]
        else:
            kappa_input[index] = TC[1]
            volume_fraction += 1
    print("volume_fraction = ", volume_fraction*1.0/dim**3)

    materials = MeshFunction("double", mesh, mesh.geometry().dim())
    for (i, cell) in enumerate(cells(mesh)):
        materials[i] = kappa_input[i / int(cell_type["num_per_cell"]) / nscale**3]
    # File("solution/materials_"+str(idx)+".pvd") << materials

    cppcode = """
    class K : public Expression
    {
    public:
        K() : Expression(1) {}
        void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
        {
            values[0] = (*materials)[cell.index];
        }
        std::shared_ptr<MeshFunction<double>> materials;
    };
    """
    kappa = Expression(cppcode=cppcode, degree=0)
    kappa.materials = materials
    return kappa[0]

log_file = open("./result.log", "w")
dump_file = open("./result.txt", "w")

for idx in range(0, 6):
    time_start = time.time()
    # Define variational problem
    kappa = GetKappa(idx)
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(0.0)
    a = dot(kappa*grad(u), grad(v))*dx
    L = f*v*dx
    # Compute solution
    u = Function(V)
    prm = GetPrm()
    solve(a == L, u, bcs, solver_parameters=prm)

    # Evaluate integral of normal gradient over top boundary
    n = FacetNormal(mesh)
    total_flux = assemble(kappa*dot(grad(u), n)*ds)
    log_file.write("total_flux = %f\n" % total_flux)
    back       = assemble(kappa*dot(grad(u), n)*ds(1))
    log_file.write("back       = %f\n" % back)
    front      = assemble(kappa*dot(grad(u), n)*ds(2))
    log_file.write("front      = %f\n" % front)
    left       = assemble(kappa*dot(grad(u), n)*ds(3))
    log_file.write("left       = %f\n" % left)
    right      = assemble(kappa*dot(grad(u), n)*ds(4))
    log_file.write("right      = %f\n" % right)
    bottom     = assemble(kappa*dot(grad(u), n)*ds(5))
    log_file.write("bottom     = %f\n" % bottom)
    top        = assemble(kappa*dot(grad(u), n)*ds(6))
    log_file.write("top        = %f\n" % top)

    # Save solution to file in VTK format
    # File("solution/solution_"+str(int(idx))+".pvd") << u

    dump_file.write("%4d%10.6f\n" % (idx, (top - bottom) / 2))
    dump_file.close()
    dump_file = open("./result.txt", "a")

    time_end = time.time()
    log_file.write("%4d time cost: %f\n\n" % (idx, time_end - time_start))
    log_file.close()
    log_file = open("./result.log", "a")

log_file.close()
dump_file.close()
