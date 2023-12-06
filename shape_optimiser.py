from firedrake import *
from firedrake.adjoint import *

continue_annotation()

# Mesh of rectangle with circle in the middle. Circle tagged with 11, Horizontal edged of rectangle
# tagged with 9 and vertical with 10
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
phi = Function(V, name="phi").interpolate(x)

# initial starting function which we later optimise on
p0 = Function(T)
p0.interpolate(1*sin(x[0]))

# for change of variables
J = grad(phi)
Jit = inv(J.T)

# setup equations to solve for u
a = (inner(u, v) + inner(Jit * grad(v), Jit * grad(u))) * det(J) * dx
L = inner(Jit("+") * normal, p0("+") * v("+")) * dS(11)
bcs = [DirichletBC(V, Constant((0,0)), "on_boundary")]
u0 = Function(V, name="u0")

# Functional to minimise -> need better one
functional = 1./2*inner(Jit * grad(w), Jit * grad(w)) * det(J) * dx

# Solve for u save into u0, here we append a cost (of u0) to our function 
dt = 1/10
for i in range(10):
    solve(a == L, u0, bcs=bcs)
    functional += Constant(dt) * (inner(u0, u0) + inner(Jit * grad(u0), Jit * grad(u0))) * det(J) * dx
    phi += dt * u0

# boundary conditions on edges of rectangle
c = 1/2 * I + 1/2
bcs = [DirichletBC(T, 1, [9]), DirichletBC(T, 0, [10])]

# Not sure what this is?
psi = TestFunction(T)
a = inner(Jit * grad(psi), Jit * grad(w)) * c * det(J) * dx
solve(a == 0, w, bcs=bcs)

# Setup for plugging into adjoint optmizer to reduce functional over p0
J = assemble(functional)
m = Control(p0)
Jhat = ReducedFunctional(J, m)

# optimise
get_working_tape().progress_bar = ProgressBar
p0_opt = minimize(Jhat)

# save outputs
out = File("minp0.pvd")
out.write(w, I, phi)