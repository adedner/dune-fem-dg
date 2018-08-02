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
# from euler import sod as problem
from euler import radialSod1Large as problem

Model, initial, x0, x1, N, name = problem(dim,gamma)

parameter.append("parameter")
parameter.append({"fem.verboserank": -1})

grid = structuredGrid(x0,x1,N)
# grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
dimR     = grid.dimension + 2
t        = 0
endTime  = 0.25
dt       = 1e-3
count    = 0
saveStep = 0.01
saveTime = saveStep

def initialize(space):
    if space.order == 0:
        return space.interpolate(initial, name='u_h')
    else:
        lagOrder = 1 # space.order
        spacelag = create.space("lagrange", space.grid, order=lagOrder, dimrange=space.dimRange)
        u_h = spacelag.interpolate(initial, name='tmp')
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

def useODESolver(polOrder=2, limiter='default'):
    global count, t, dt, saveTime
    polOrder = polOrder
    if False:
        # needs change in dune/fem-dg/operator/dg/passtraits.hh
        space = create.space("dglegendre", grid, order=polOrder, dimrange=dimR, hierarchical=False)
    else:
        space = create.space("dgonb", grid, order=polOrder, dimrange=dimR)
    u_h = initialize(space)
    rho, v, p = Model.toPrim(u_h)
    operator = createFemDGSolver( Model, space, limiter=limiter )
    # operator.setTimeStepSize(dt)

    start = time.time()
    # operator.applyLimiter( u_h );
    print(grid.size(0))
    grid.writeVTK(name,
        pointdata=[u_h],
        # celldata={"density":rho, "pressure":p}, # bug: density not shown correctly
        celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
        cellvector={"velocity":v},
        number=count, subsampling=2)
    tcount = 0
    while t < endTime:
        operator.solve(target=u_h)
        # operator.applyLimiter( u_h );
        dt = operator.deltaT()
        t += dt
        tcount += 1
        if t > saveTime:
            print('dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
            count += 1
            # rho, v, p = Model.toPrim(u_h) # is this needed - works for me # without and it should...
            grid.writeVTK(name,
                pointdata=[u_h],
                celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
                cellvector={"velocity":v},
                number=count, subsampling=2)
            saveTime += saveStep
    grid.writeVTK(name,
        pointdata=[u_h],
        celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
        cellvector={"velocity":v},
        number=count, subsampling=2)
    print("time loop:",time.time()-start)
    print("number of time steps ", tcount)


if True:
    grid = structuredGrid(x0,x1,N)
    # grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    useODESolver(2,'default')      # third order with limiter
elif False:
    N = [n*10 for n in N]
    grid = structuredGrid(x0,x1,N)
    # grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    useODESolver(0,None)           # FV scheme
elif False:
    N = [n*10 for n in N]
    grid = structuredGrid(x0,x1,N)
    # grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    useODESolver(0,'default')      # FV scheme with limiter
