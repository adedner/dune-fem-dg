import time
from dune.grid import structuredGrid, cartesianDomain
from dune.fem import parameter
import dune.create as create
from dune.femdg import createFemDGSolver
# from ufl import *

from shallowwater import leVeque as problem

Model, initial, x0,x1,N, endTime, name = problem()

parameter.append({"fem.verboserank": -1})
parameter.append("parameter")

grid = structuredGrid(x0,x1,N)
dimR     = grid.dimension + 1
t        = 0
count    = 0
saveStep = 0.01
saveTime = saveStep

def initialize(space):
    return space.interpolate(Model.toCons(initial), name='u_h')

def useODESolver(polOrder=2, limiter='default'):
    global count, t, saveTime
    polOrder = polOrder
    space = create.space("dgonb", grid, order=polOrder, dimrange=dimR)
    u_h = initialize(space)
    eta, v = Model.toPrim(u_h)
    operator = createFemDGSolver( Model, space, limiter=limiter )

    operator.applyLimiter( u_h );
    print("number of elements: ",grid.size(0),flush=True)
    grid.writeVTK(name,
        celldata=[u_h],
        pointdata={"freeSurface":eta},
        cellvector={"velocity":v},
        number=count, subsampling=2)
    start = time.time()
    tcount = 0
    while t < endTime:
        operator.solve(target=u_h)
        dt = operator.deltaT()
        t += dt
        tcount += 1
        if tcount%100 == 0:
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
        if t > saveTime:
            count += 1
            grid.writeVTK(name,
                celldata=[u_h],
                pointdata={"freeSurface":eta},
                cellvector={"velocity":v},
                number=count, subsampling=2)
            saveTime += saveStep
    print("time loop:",time.time()-start)
    print("number of time steps ", tcount)
    grid.writeVTK(name,
        celldata=[u_h],
        pointdata={"freeSurface":eta},
        cellvector={"velocity":v},
        number=count, subsampling=2)

useODESolver(2,'default')      # third order with limiter
