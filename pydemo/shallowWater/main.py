import time
from dune.grid import structuredGrid, cartesianDomain
from dune.fem import parameter
import dune.create as create
from dune.femdg import createFemDGSolver

from shallowwater import leVeque as problem

Model, initial, x0,x1,N, endTime, name = problem()

parameter.append({"fem.verboserank": -1})
parameter.append("parameter")

grid     = structuredGrid(x0,x1,N)
dimR     = Model.dimension
t        = 0
count    = 0
saveStep = 0.01
saveTime = saveStep
polOrder = 2
limiter  = "default"

space = create.space("dgonb", grid, order=polOrder, dimrange=dimR)
u_h = space.interpolate(Model.toCons(initial), name='u_h')
eta, v = Model.toPrim(u_h)
operator = createFemDGSolver( Model, space, limiter=limiter )

operator.applyLimiter( u_h );
print("number of elements: ",grid.size(0),flush=True)
vtk = grid.writeVTK(name, subsampling=2, write=False,
           celldata=[u_h],
           pointdata={"freeSurface":eta},
           cellvector={"velocity":v}
        )
vtk.write(name, count)
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
        vtk.write(name, count)
        saveTime += saveStep
print("time loop:",time.time()-start)
print("number of time steps ", tcount)
vtk.write(name, count)
