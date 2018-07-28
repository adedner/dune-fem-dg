import time
from dune.grid import structuredGrid, cartesianDomain
from dune.fem import parameter
import dune.create as create
from dune.models.elliptic.formfiles import loadModels
from dune.femdg import createFemDGSolver

# from scalar import riemanProblem as problem
# from scalar import constantTransport as probelm
# from scalar import cosProblem as problem
from scalar import cosProblem as problem
# from scalar import burgersShock as problem
# from scalar import burgersVW as problem
# from scalar import burgersStationaryShock as problem

Model, initial = problem()


parameter.append("parameter")
parameter.append({"fem.verboserank": -1})

x0,x1,N = [-1, 0], [1, 0.1], [50, 7]
grid = structuredGrid(x0,x1,N)
# grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
dimR     = 1
t        = 0
saveStep = 0.01
saveTime = saveStep
count    = 0
endTime  = 1.

def useODESolver():
    global count, t, saveTime

    spaceName = "dgonb"
    polOrder = 2
    space = create.space(spaceName, grid, order=polOrder, dimrange=dimR)
    u_h   = space.interpolate(initial, name='u_h')
    u_h_n = u_h.copy(name="previous")
    operator = createFemDGSolver( Model, space )

    start = time.time()
    grid.writeVTK(Model.name, celldata=[u_h], number=count)
    while t < endTime:
        u_h_n.assign(u_h)
        operator(u_h_n, u_h)
        dt = operator.deltaT()
        t += dt
        if t > saveTime:
            print('dt = ', dt, 'time = ',t, 'count = ',count )
            count += 1
            grid.writeVTK(Model.name, celldata=[u_h], number=count)
            saveTime += saveStep
    grid.writeVTK(Model.name, celldata=[u_h], number=count)
    print("time loop:",time.time()-start)

useODESolver()
