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
functional = (1-I)*1./2*inner(Jit * grad(w), Jit * grad(w)) * det(J) * dx
functional = inner(x,x)*dx

# Solve for u save into u0, here we append a cost (of u0) to our function 
dt = 1/10
for i in range(10):
    solve(a == L, u0, bcs=bcs)
    functional += Constant(dt/100) * (inner(u0, u0) + inner(Jit * grad(u0), Jit * grad(u0))) * det(J) * dx
    dphi += dt * u0

# boundary conditions on edges of rectangle
bcs = [DirichletBC(T, 1, [9]), DirichletBC(T, 0, [10])]

# shape dependent equation
psi = TestFunction(T)
c = 100 * I + 1
a = (psi*w + inner(Jit * grad(psi), Jit * grad(w))) * c * det(J) * dx
solve(a == 0, w, bcs=bcs)

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
