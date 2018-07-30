import time
from dune.grid import structuredGrid, cartesianDomain
from dune.fem import parameter
import dune.create as create
from dune.models.elliptic.formfiles import loadModels
from llf import NumFlux
from dune.femdg import createFemDGSolver
from ufl import *

gamma = 1.4
dim = 2
from euler import sod as problem
# from euler import radialSod3 as problem

Model, initial, x0, x1, N, name = problem(dim,gamma)

parameter.append("parameter")
parameter.append({"fem.verboserank": -1})

grid = structuredGrid(x0,x1,N)
# grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
dimR      = 4
t = 0
dt = 1e-5
saveStep = 0.01
saveTime = saveStep
count = 0
endTime = 0.4

def initialize(space):
    lagOrder = 1 # space.order
    lag = create.space("lagrange", space.grid, order=1, dimrange=space.dimRange)
    u_h = lag.interpolate(initial, name='tmp')
    return space.interpolate(u_h, name='u_h')

def useGalerkinOp():
    global count, t, dt, saveTime
    # a full fv implementation using UFL and the galerkin operator
    # some sign error or something else leading to wrong results
    # still needs fixing
    spaceName = "finitevolume"
    space = create.space(spaceName, grid, dimrange=dimR)

    n = FacetNormal(space.cell())

    u_h   = space.interpolate(initial, name='u_h')
    u_h_n = u_h.copy(name="previous")

    fullModel = inner( Model.F_c(u), grad(v) ) * dx -\
                inner( NumFlux(Model,u,n,dt), v('+')-v('-')) * dS -\
                inner( Model.F_c(u)*n, v) * ds
    operator = create.operator("galerkin",
                 inner(u-u_h_n,v)*dx == dt*replace(fullModel,{u:u_h_n}),
                 space )

    operator.model.dt = 1e-5

    start = time.time()
    grid.writeVTK(name, pointdata=[u_h], number=count)
    while t < endTime:
        u_h_n.assign(u_h)
        operator.solve(target=u_h)
        t += operator.model.dt
        if t > saveTime:
            count += 1
            grid.writeVTK(name, pointdata=[u_h], number=count)
            saveTime += saveStep
    grid.writeVTK(name, pointdata=[u_h], number=count)
    print("time loop:",time.time()-start)

def useODESolver():
    global count, t, dt, saveTime
    spaceName = "dgonb"
    polOrder = 2
    space = create.space(spaceName, grid, order=polOrder, dimrange=dimR)
    u_h = initialize(space)
    rho, v, p = Model.toPrim(u_h)
    operator = createFemDGSolver( Model, space, limiter=None )
    # operator.setTimeStepSize(dt)

    start = time.time()
    grid.writeVTK(name,
        pointdata=[u_h],
        celldata={"density":rho, "pressure":p},
        cellvector={"velocity":v},
        number=count)
    while t < endTime:
        # operator.applyLimiter( u_h );
        operator.solve(target=u_h)
        dt = operator.deltaT()
        t += dt
        if t > saveTime:
            print('dt = ', dt, 'time = ',t, 'count = ',count )
            count += 1
            grid.writeVTK(name,
                pointdata=[u_h],
                celldata={"density":rho, "pressure":p},
                cellvector={"velocity":v},
                number=count)
            saveTime += saveStep
    grid.writeVTK(name,
        pointdata=[u_h],
        celldata={"density":rho, "pressure":p},
        cellvector={"velocity":v},
        number=count)
    print("time loop:",time.time()-start)

useODESolver()
