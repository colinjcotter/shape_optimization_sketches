from firedrake import *
from firedrake.adjoint import *

continue_annotation()

# Mesh of rectangle with circle in the middle. Circle tagged with 11, Horizontal edged of rectangle
# tagged with 9 and vertical with 10
# gmsh -2 1.geo -format msh2
mesh = Mesh('1.msh')

S = FunctionSpace(mesh, "DG", 0)
V = VectorFunctionSpace(mesh, "CG", 1)
T = FunctionSpace(mesh, "CG", 1)

x = SpatialCoordinate(mesh)
I = Function(S, name="indicator")
I.interpolate(conditional(sqrt((x[0]-5)**2 + (x[1]+4)**2) < 2, 1, 0))

n = FacetNormal(mesh)
normal = I("+")*n("+") + I("-")*n("-")

v = TestFunction(V)
u = TrialFunction(V)
w = Function(T, name="w")

# reparametrisation to not have to move the mesh
dphi = Function(V, name="dphi")
phi = x + dphi

# initial starting function which we later optimise on
p0 = Function(T)

# for change of variables
J = grad(phi)
Jit = inv(J.T)

# setup equations to solve for u
a = (inner(u, v) + inner(Jit * grad(v), Jit * grad(u))) * det(J) * dx
L = inner(Jit("+") * normal, p0("+") * v("+")) * dS(11)
bcs = [DirichletBC(V, Constant((0,0)), "on_boundary")]
u0 = Function(V, name="u0")

# Functional to minimise -> need better one
f = (x[0]+dphi[0]-0.5)**2 + (x[1]+dphi[1]+4)**2 - 0.7
functional = f*I*det(J)*dx

# Solve for u save into u0, here we append a cost (of u0) to our function 
dt = 1/10
for i in range(10):
    solve(a == L, u0, bcs=bcs)
    dphi += dt * u0

# Setup for plugging into adjoint optmizer to reduce functional over p0
J = assemble(functional)
m = Control(p0)
Jhat = ReducedFunctional(J, m)

# optimise
get_working_tape().progress_bar = ProgressBar
p0_opt = minimize(Jhat)
print(p0_opt.dat.data)

# save outputs
out = File("minp0.pvd")
out.write(w, I, dphi)
