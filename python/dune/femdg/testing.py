import time, math
from dune.grid import structuredGrid
import dune.create as create
from dune.fem.function import integrate
from dune.femdg import createFemDGSolver
from ufl import dot

def run(Model, initial, x0,x1,N, endTime, name, exact,
        polOrder, limiter="default", startLevel=0,
        primitive=None, saveStep=None, subsamp=0):
    grid     = structuredGrid(x0,x1,N)
    grid.hierarchicalGrid.globalRefine(startLevel)
    dimR     = Model.dimension
    t        = 0
    count    = 0
    saveTime = saveStep

    space = create.space("dgonb", grid, order=polOrder, dimrange=dimR)
    u_h = space.interpolate(Model.toCons(initial), name='u_h')
    operator = createFemDGSolver( Model, space, limiter=limiter )

    operator.applyLimiter( u_h );
    print("number of elements: ",grid.size(0),flush=True)
    if saveStep is not None:
        vtk = grid.writeVTK(name, subsampling=subsamp, write=False,
               celldata=[u_h],
               pointdata=primitive(Model,u_h) if primitive else None,
               cellvector={"velocity":Model.velo(u_h)}
            )
        vtk.write(name, count)
    start = time.time()
    tcount = 0
    while t < endTime:
        operator.solve(target=u_h)
        if math.isnan( u_h.scalarProductDofs( u_h ) ):
            print('ERROR: dofs invalid t =', t)
            exit(0)
        dt = operator.deltaT()
        t += dt
        tcount += 1
        if True: # tcount%100 == 0:
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
        if saveStep is not None and t > saveTime:
            count += 1
            vtk.write(name, count)
            saveTime += saveStep
    print("time loop:",time.time()-start)
    print("number of time steps ", tcount)
    vtk.write(name, count)

    if exact is not None:
        error = integrate( grid, dot(u_h-exact(t),u_h-exact(t)), order=5 )
        print("error:", math.sqrt(error) )
