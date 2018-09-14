import time, math
from dune.grid import structuredGrid, cartesianDomain, OutputType
import dune.create as create
from dune.fem.function import integrate
from dune.ufl import NamedConstant
from ufl import dot, SpatialCoordinate

def run(Model, initial, x0,x1,N, endTime, name, exact,
        polOrder, limiter="default", startLevel=0,
        primitive=None, saveStep=None, subsamp=0,
        dt=None,grid="yasp",threading=True,
        parameters={}):
    periodic=[True,]*len(x0)
    if hasattr(Model,"boundary"):
        bnd=set()
        for b in Model.boundary:
            bnd.update(b)
        for i in range(len(x0)):
            if 2*i+1 in bnd:
                assert(2*i+2 in bnd)
                periodic[i] = False
    if any(bnd):
        print("Setting periodic boundaries",periodic,flush=True)

    domain   = cartesianDomain(x0,x1,N,periodic=periodic,overlap=0)
    grid     = create.grid(grid,domain)
    grid.hierarchicalGrid.globalRefine(startLevel)
    dimR     = Model.dimension
    t        = 0
    count    = 0
    saveTime = saveStep

    space = create.space("dgonb", grid, order=polOrder, dimrange=dimR)
    u_h = space.interpolate(initial, name='u_h')
    operator = create.scheme("femDG",Model, space, limiter=limiter, threading=True, parameters=parameters )
    operator.applyLimiter( u_h );
    print("number of elements: ",grid.size(0),flush=True)
    if saveStep is not None:
        x = SpatialCoordinate(space.cell())
        tc = NamedConstant(space.cell(),"time")
        try:
            velo = [create.function("ufl",space.grid, ufl=Model.velocity(tc,x,u_h), order=2, name="velocity")]
        except AttributeError:
            velo = None
        vtk = grid.writeVTK(name, subsampling=subsamp, write=False,
               celldata=[u_h],
               pointdata=primitive(Model,u_h) if primitive else None,
               cellvector=velo
            )
        try:
            velo[0].setConstant("time",[t])
        except:
            pass
        vtk.write(name, count) # , outputType=OutputType.appendedraw)
    start = time.time()
    tcount = 0
    while t < endTime:
        if dt is not None:
            operator.setTimeStepSize(dt)
        operator.step(target=u_h)
        dt = operator.deltaT()
        if math.isnan( u_h.scalarProductDofs( u_h ) ):
            grid.writeVTK(name, subsampling=subsamp, celldata=[u_h])
            print('ERROR: dofs invalid t =', t,flush=True)
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
            exit(0)
        t += dt
        tcount += 1
        if tcount%100 == 0:
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
        if saveStep is not None and t > saveTime:
            count += 1
            try:
                velo[0].setConstant("time",[t])
            except:
                pass
            vtk.write(name, count, outputType=OutputType.appendedraw)
            saveTime += saveStep

    print("time loop:",time.time()-start,flush=True)
    print("number of time steps ", tcount,flush=True)
    if saveStep is not None:
        try:
            velo[0].setConstant("time",[t])
        except:
            pass
        vtk.write(name, count, outputType=OutputType.appendedraw)

    # output the final result and compute error (if exact is available)
    if exact is not None:
        grid.writeVTK(name, subsampling=subsamp,
                celldata={"solution":u_h, "exact":exact(t)})
        error = integrate( grid, dot(u_h-exact(t),u_h-exact(t)), order=5 )
        print("error:", math.sqrt(error),flush=True )
    else:
        grid.writeVTK(name, subsampling=subsamp, celldata=[u_h])
