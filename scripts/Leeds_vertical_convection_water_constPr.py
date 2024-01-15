# Leeds_vertical_convection_vortices.py
# incompressible fluid convection, steady-state, hot on left upright, cold on right upright
# this is based on example from https://www.firedrakeproject.org/demos/rayleigh-benard.py
# with small changes by Ed Threlfall, January 2023

# converted to allow Prandtl number to depend on temperature
# Ed Threlfall, December 2023
# QoI is Nusselt number = heat flux left-right, called Nusselt in this script
# de-dimensionalization done by rescaling time to inverse of free-fall velocity

# this example enables examination of T-dependent Pr for water but is currently set up with constant Pr of 7.0

from firedrake import *

from firedrake.petsc import PETSc
import math
import numpy as np

height = 10
meshres = 64

#M = Mesh("refined_unit_square.msh")  # can be used with boundary-refined mesh for larger Ra
M = RectangleMesh(meshres, meshres*height, 1, height, quadrilateral=True)

# choice of finite element space is Taylor-Hood elements:
V = VectorFunctionSpace(M, "CG", 2)
W = FunctionSpace(M, "CG", 1)
Q = FunctionSpace(M, "CG", 1)
Z = V * W * Q

upT = Function(Z)
u, p, T = split(upT)
v, q, S = TestFunctions(Z)

x, y = SpatialCoordinate(M)

Ra=Constant(1.0e4)  # Rayleigh number

amp_0 = 7.0  # Prandtl number for water (room temp) - gets overwritten later in this example
amp_coeffs = np.zeros(3)
amp_coeffs[0] =  amp_0

# coefficients for order-2 Taylor series, Pr(T)
amp_coeffs[1] = 0.0  # for water try -5.0 for initial attempt, means Pr is 2.0 at hot edge and 7.0 at cold edge (room temp)
amp_coeffs[2] = 0.0  # use 0.0 for initial attempt
# end of parameters to vary

g = Constant((0, 1))  # buoyancy force

bcs = [
    DirichletBC(Z.sub(0), Constant((0, 0)), (1, 2, 3, 4)),
    DirichletBC(Z.sub(2), Constant(1.0), (1,)),
    DirichletBC(Z.sub(2), Constant(0.0), (2,)),
]

# Like Navier-Stokes, the pressure is only defined up to a constant.::
nullspace = MixedVectorSpaceBasis(
    Z, [Z.sub(0), VectorSpaceBasis(constant=True), Z.sub(2)])

from firedrake.petsc import PETSc

# solve
for i in range (0,1):  # continuation optional, no continuation means run once (use range (0,1))

   #amp_coeffs[0] = pow(10,-i*0.1)

   F = (
       sqrt((amp_coeffs[0]+amp_coeffs[1]*T+amp_coeffs[2]*T*T)/Ra) * inner(grad(u), grad(v))*dx
       + inner(dot(grad(u), u), v)*dx
       - inner(p, div(v))*dx
       - inner(T*g, v)*dx
       + inner(div(u), q)*dx
       + inner(dot(grad(T), u), S)*dx
       + (1/sqrt((amp_coeffs[0]+amp_coeffs[1]*T+amp_coeffs[2]*T*T)*Ra))*inner(grad(T), grad(S))*dx
   )

   # Nusselt is ratio of heat flux with flow to heat flux with no flow, hence normalization:

   T_n = Function(Q)
   S_n = TestFunction(Q)

   bcs_n = [ DirichletBC(Q, Constant(1.0), (1,)), DirichletBC(Q, Constant(0.0), (2,))]

   F_n = (
        (1/sqrt((amp_coeffs[0]+amp_coeffs[1]*T_n+amp_coeffs[2]*T_n*T_n)*Ra))*inner(grad(T_n), grad(S_n))*dx
   )

   print("starting run ... \n")
   print("i="+str(int(i))+"\n")
   print("Ra="+str(float(Ra))+"\n")
   print("amp_coefficients[0]="+str(float(amp_coeffs[0]))+"\n")
   try:
      solve(F == 0, upT, bcs=bcs, nullspace=nullspace,
            solver_parameters={"mat_type": "aij",
                               "snes_monitor": None,
                               "ksp_type": "gmres",
                               "pc_type": "lu",
                               "pc_factor_mat_solver_type": "mumps"})

      solve(F_n == 0, T_n, bcs=bcs_n,
            solver_parameters={"mat_type": "aij",
                               "snes_monitor": None,
                               "ksp_type": "gmres",
                               "pc_type": "lu",
                               "pc_factor_mat_solver_type": "mumps"})

   except:
      print("something bad happened on run!")
      quit()  # optionally quit
      #upT.sub(2).interpolate(Constant(0.0))  # this is done to make Nu appear zero for runs that have failed

# evaluate QoI = Nusselt number
normL = Function(V)
normL = Constant((-1.0,0.0))
fluxL = assemble(inner(normL, grad(T))*ds(1))
normR = Function(V)
normR = Constant((1.0,0.0))
fluxR = assemble(inner(normR, grad(T))*ds(2))
PrL = amp_0+amp_coeffs[1]*1+amp_coeffs[2]*1*1

fluxL_normalize = assemble(inner(normL, grad(T_n))*ds(1))
fluxR_normalize = assemble(inner(normR, grad(T_n))*ds(2))

PrR = amp_0
Nusselt = (fluxL/sqrt(PrL)-fluxR/sqrt(PrR))/(fluxL_normalize/sqrt(PrL)-fluxR_normalize/sqrt(PrR))
print("finished run with log_10 Ra "+str(float(math.log10(Ra)))+" and Nusselt number "+str(float(Nusselt)))

# these ought to be identical for steady-state solutions
print("heat flux left "+str(float(fluxL/sqrt(PrL))))
print("heat flux rght "+str(float(fluxR/sqrt(PrR))))

# as should these
print("normalization heat flux left "+str(float(fluxL_normalize/sqrt(PrL))))
print("normalization heat flux rght "+str(float(fluxR_normalize/sqrt(PrR))))

# generate ParaView output, comment out if desired
u, p, T = upT.split()
u.rename("Velocity")
p.rename("Pressure")
T.rename("Temperature")
Pr = Function(Z.sub(1))  # Prandtl number: now a function of T (temperature)
Pr.interpolate((amp_coeffs[0]+amp_coeffs[1]*T+amp_coeffs[2]*T*T))
Pr.rename("Prandtl number")
File("Leeds_vertical_convection_water_constPr.pvd").write(u,p,T,Pr)

Pr_n = Function(Z.sub(1))  # Prandtl number: now a function of T (temperature)
Pr_n.interpolate((amp_coeffs[0]+amp_coeffs[1]*T_n+amp_coeffs[2]*T_n*T_n))
File("Leeds_vertical_convection_water__constPr_norm.pvd").write(T_n,Pr_n)







