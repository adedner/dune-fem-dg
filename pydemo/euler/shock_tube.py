import time
from dune.grid import structuredGrid, cartesianDomain
from dune.fem import parameter
import dune.create as create
from dune.models.elliptic.formfiles import loadModels
from llf import NumFlux
from dune.femdg import createFemDGSolver
from ufl import *

from euler import Model
from euler import Sod as Initial

parameter.append("parameter")
parameter.append({"fem.verboserank": -1})

x0,x1,N = [-1, 0], [1, 0.1], [20, 5]
grid = structuredGrid(x0,x1,N)
# grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
dimR      = 4
t = 0
dt = 1e-5
saveStep = 0.01
saveTime = saveStep
count = 0
endTime = 0.4

def useGalerkinOp():
    global count, t, dt, saveTime
    # a full fv implementation using UFL and the galerkin operator
    # some sign error or something else leading to wrong results
    # still needs fixing
    spaceName = "finitevolume"
    space = create.space(spaceName, grid, dimrange=dimR)

    n = FacetNormal(space.cell())

    u_h   = space.interpolate(Initial(), name='u_h')
    u_h_n = u_h.copy(name="previous")

    fullModel = inner( Model.F_c(u), grad(v) ) * dx -\
                inner( NumFlux(Model,u,n,dt), v('+')-v('-')) * dS -\
                inner( Model.F_c(u)*n, v) * ds
    operator = create.operator("galerkin",
                 inner(u-u_h_n,v)*dx == dt*replace(fullModel,{u:u_h_n}),
                 space )

    operator.model.dt = 1e-5

    start = time.time()
    grid.writeVTK('sod', pointdata=[u_h], number=count)
    while t < 0.4:
        u_h_n.assign(u_h)
        operator.solve(target=u_h)
        t += operator.model.dt
        if t > saveTime:
            count += 1
            grid.writeVTK('sod', pointdata=[u_h], number=count)
            saveTime += saveStep
    grid.writeVTK('sod', pointdata=[u_h], number=count)
    print("time loop:",time.time()-start)

def useODESolver():
    global count, t, dt, saveTime
    spaceName = "dgonb"
    polOrder = 2
    space = create.space(spaceName, grid, order=polOrder, dimrange=dimR)
    u_h   = space.interpolate(Initial(), name='u_h')
    u_h_n = u_h.copy(name="previous")
    operator = createFemDGSolver( Model, space )
    # operator.setTimeStepSize(dt)

    start = time.time()
    grid.writeVTK('sod', pointdata=[u_h], number=count)
    while t < 0.4:
        u_h_n.assign(u_h)
        operator(u_h_n, u_h)
        dt = operator.deltaT()
        t += dt
        if t > saveTime:
            print('dt = ', dt, 'time = ',t, 'count = ',count )
            count += 1
            grid.writeVTK('sod', pointdata=[u_h], number=count)
            saveTime += saveStep
    grid.writeVTK('sod', pointdata=[u_h], number=count)
    print("time loop:",time.time()-start)

useODESolver()
