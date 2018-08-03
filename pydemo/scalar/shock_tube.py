import time, math
from dune.grid import structuredGrid, cartesianDomain
from dune.fem import parameter
from dune.fem.function import integrate
import dune.create as create
from dune.models.elliptic.formfiles import loadModels
from dune.femdg import createFemDGSolver
from ufl import dot

# from scalar import shockTransport as problem
# from scalar import constantTransport as probelm
# from scalar import sinProblem as problem
from scalar import pulse as problem

Model, initial, x0,x1,N, endTime, name, exact = problem()

parameter.append({"fem.verboserank": 0})
parameter.append("parameter")

grid = structuredGrid(x0,x1,N)

# grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
dimR     = 1
t        = 0
dt       = 0.005
saveStep = 0.05
saveTime = saveStep
count    = 0

def useODESolver():
    global count, t, dt, saveTime

    spaceName = "dgonb"
    polOrder = 2
    space = create.space(spaceName, grid, order=polOrder, dimrange=dimR)
    u_h   = space.interpolate(initial, name='u_h')
    operator = createFemDGSolver( Model, space, limiter=None, diffusionScheme="cdg2" )

    operator.setTimeStepSize(dt)

    start = time.time()
    grid.writeVTK(name, celldata=[u_h], number=count, subsampling=2)
    while t < endTime:
        operator.solve(target=u_h)
        dt = operator.deltaT()
        t += dt
        if t > saveTime:
            print('dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
            count += 1
            grid.writeVTK(name, celldata=[u_h], number=count, subsampling=2)
            saveTime += saveStep
    grid.writeVTK(name, celldata=[u_h], number=count, subsampling=2)
    print("time loop:",time.time()-start)
    if exact is not None:
        error = integrate( grid, dot(u_h-exact(t),u_h-exact(t)), order=5 )
        print("error:", math.sqrt(error) )

useODESolver()
