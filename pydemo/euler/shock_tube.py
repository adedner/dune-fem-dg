import time
from dune.grid import structuredGrid
from dune.fem import parameter
import dune.create as create
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, exp,\
                FacetNormal, dx, grad, inner, as_vector, replace, sqrt,\
                conditional, as_matrix, Max, dS, ds, avg, jump, outer
from dune.ufl import NamedConstant
from dune.femdg import createFemDGSolver

parameter.append({"fem.verboserank": -1})

grid = structuredGrid([-1, 0], [1, 0.1], [200, 10])

order     = 0
dimR      = 4
spaceName = "dgonb"
space = create.space(spaceName, grid, dimrange=dimR, order=order)

import euler

Euler = euler.Model
NumFLux = euler.NumFlux

u = TrialFunction(space)
v = TestFunction(space)
dt = NamedConstant(triangle, "dt")    # time step
x = SpatialCoordinate(space.cell())
n = FacetNormal(space.cell())

UL = as_vector( Euler.toCons([1,0,0,1]) )
UR = as_vector( Euler.toCons([0.125,0,0,0.1]) )
initial = conditional(x[0]<0,UL,UR)

u_h = space.interpolate(initial, name='u_h')
u_h_n = u_h.copy(name="previous")

# here a full implementation using UFL and the galerkin operator
fullModel = inner(u_h_n,v)*dx-euler.F +\
            dt*inner(LLFFlux.flux(Euler,u,n,dt),avg(v))*dS +\
            dt*inner(Euler.F_c(u)*n,v)*ds
# operator = create.operator("galerkin", fullModel, space )
operator = createFemDGSolver( "euler", space, fullModel )

operator.model.dt = 5e-4

t = 0
saveStep = 0.01
saveTime = saveStep
count = 0
grid.writeVTK('sod', pointdata=[u_h], number=count)

start = time.time()
while t < 0.4:
    u_h_n.assign(u_h)
    operator(u_h_n,u_h)
    t += operator.model.dt
    if t > saveTime:
        count += 1
        grid.writeVTK('sod', pointdata=[u_h], number=count)
        saveTime += saveStep
grid.writeVTK('sod', pointdata=[u_h], number=count)
print("time loop:",time.time()-start)
