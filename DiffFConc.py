"""
Fenics program on an attempt to create a way to model concentration dependent
diffusivity. This will be a NON-LINEAR solver
"""

from __future__ import print_function
from fenics import *
import numpy as np
from subprocess import call

T = 5.0
num_steps = 10
dt = T/num_steps

# Define the non-linear coefficient Diff
def Diff(u):
	"Return nonlinear coefficient"
	#return 5.0E-6 * (1.0 - 0.36*u + ((0.36*u)**2.0)/2.0)
	return 5.0E-6 * np.exp(-0.36*u) #np.exp(-0.36*u)



# Create mesh and define function space
nx = 96
ny = 96
mesh = RectangleMesh(Point(-1,0), Point(1,2),nx,ny)
V = FunctionSpace(mesh, 'P', 1)

# Define the tolerance of the boundary condition
tol = 1.0E-14

# Define the boundary condition
# Including mixed Neumann (0 at boundaries) and Dirichlet conditions
u_D = Expression('3.0 * (1 - exp(-5000 * t))', degree=2, t=0) # POSSIBLE ISSUE IN DEGREE

def boundary_D(x, on_boundary):
	if on_boundary:
		if near(x[1], 0, tol) and near(x[0],0,0.8):
			return True
		else:
			return False
	else:
		return False


bc = DirichletBC(V, u_D, boundary_D)

# Define the initial value
u_n = interpolate(u_D, V)


# Define the variational problem
u = Function(V) # Note: NOT TRIALFUNCTION!
v = TestFunction(V)
f = Constant(0.0)

#a = (u*v + dt*Diff(u)*dot(grad(u), grad(v)))*dx
#L = (u_n + dt*f)*v*dx
F = (u*v + dt*Diff(u)*dot(grad(u),grad(v)))*dx - (u_n + dt*f)*v*dx

# Create VTK file for saving solution
vtkfile = File('DiffFConc/solution.pvd')

# Time stepping
u = Function(V)
t = 0.0

for n in range(num_steps):

	# Update current time
	t += dt
	u_D.t = t

	# Compute solution
	solve(F == 0, u, bc)

	# Save to file and plot solution
	vtkfile << (u,t)

	# Print time
	print('t= %2.f' % t)

	# Update previous solution
	u_n.assign(u)

