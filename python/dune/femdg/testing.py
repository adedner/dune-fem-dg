import time, math, sys
from dune.grid import structuredGrid, cartesianDomain, OutputType
import dune.create as create
from dune.fem.function import integrate
from dune.ufl import Constant
from ufl import dot, SpatialCoordinate
from dune.femdg import femDGOperator, femDGSolver, rungeKuttaSolver

from collections import namedtuple

Parameters = namedtuple("TestingParameters",
                        ["Model", "initial", "domain", "endTime", "name", "exact"])
Parameters.__new__.__defaults__ = (None,None) # defaults for name, exact

def run(Model, initial, domain, endTime, name, exact,
        polOrder, limiter="default", startLevel=0,
        primitive=None, saveStep=None, subsamp=0,
        dt=None,grid="yasp", space="dgonb", threading=True,
        parameters={}):
    print("*************************************")
    print("**** Running simulation",name)
    print("*************************************")
    try: # passed in a [xL,xR,N] tripple
        x0,x1,N = domain
        periodic=[True,]*len(x0)
        if hasattr(Model,"boundary"):
            bnd=set()
            for b in Model.boundary:
                bnd.update(b)
            for i in range(len(x0)):
                if 2*i+1 in bnd:
                    assert(2*i+2 in bnd)
                    periodic[i] = False
        print("Setting periodic boundaries",periodic,flush=True)
        # create domain and grid
        domain   = cartesianDomain(x0,x1,N,periodic=periodic,overlap=0)
        grid     = create.grid(grid,domain)
        # initial refinement of grid
        grid.hierarchicalGrid.globalRefine(startLevel)
    except TypeError: # assume the 'domain' is already a gridview
        grid = domain

    dimR     = Model.dimRange
    t        = 0
    count    = 0
    saveTime = saveStep

    # create discrete function space
    space = create.space( space, grid, order=polOrder, dimRange=dimR)
    # create and initialize solution
    u_h = space.interpolate(initial, name='u_h')
    # create solution scheme, i.e. operator and ODE solver

    # create DG operator based on Model
    operator = femDGOperator(Model, space, limiter=limiter, threading=True, parameters=parameters )
    # create Runge-Kutta solver
    rkScheme = rungeKuttaSolver( operator, parameters=parameters )

    # limit initial data if necessary
    operator.applyLimiter( u_h );

    print("number of elements: ",grid.size(0),flush=True)

    # preparation for output
    if saveStep is not None:
        x = SpatialCoordinate(space.cell())
        tc = Constant(0.0,"time")
        try:
            velo = [create.function("ufl",space.grid, ufl=Model.velocity(tc,x,u_h), order=2, name="velocity")]
        except AttributeError:
            velo = None
        if name is not None:
            vtk = grid.writeVTK(name, subsampling=subsamp, write=False,
                   celldata=[u_h],
                   pointdata=primitive(Model,u_h) if primitive else None,
                   cellvector=velo
                )
        else:
            vtk = None
        try:
            velo[0].setConstant("time",[t])
        except:
            pass
        if vtk is not None:
            vtk.write(name, count) # , outputType=OutputType.appendedraw)

    # measure CPU time
    start = time.time()

    tcount = 0
    # time loop
    # set time step size to ODE solver
    if dt is not None:
        rkScheme.setTimeStepSize(dt)
    print("Start solving")
    while t < endTime:
        # solver time step
        assert not math.isnan( u_h.scalarProductDofs( u_h ) )
        rkScheme.solve(u_h)

        # limit solution if necessary
        operator.applyLimiter( u_h )

        # obtain new time step size
        dt = rkScheme.deltaT()
        # check that solution is meaningful
        if name is not None and math.isnan( u_h.scalarProductDofs( u_h ) ):
            grid.writeVTK(name, subsampling=subsamp, celldata=[u_h])
            print('ERROR: dofs invalid t =', t,flush=True)
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
            sys.exit(1)
        # increment time and time step counter
        t += dt
        tcount += 1

        # output
        if True: # tcount%100 == 0:
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
        if saveStep is not None and t > saveTime:
            count += 1
            try:
                velo[0].setConstant("time",[t])
            except:
                pass
            if vtk is not None:
                vtk.write(name, count, outputType=OutputType.appendedraw)
            saveTime += saveStep

    print("time loop:",time.time()-start,flush=True)
    print("number of time steps ", tcount,flush=True)
    if name is not None and saveStep is not None:
        try:
            velo[0].setConstant("time",[t])
        except:
            pass
        vtk.write(name, count, outputType=OutputType.appendedraw)

    # output the final result and compute error (if exact is available)
    if exact is not None and name is not None:
        grid.writeVTK(name, subsampling=subsamp,
                celldata=[u_h], pointdata={"exact":exact(t)})
        error = integrate( grid, dot(u_h-exact(t),u_h-exact(t)), order=5 )
        print("error:", math.sqrt(error),flush=True )
    elif name is not None:
        grid.writeVTK(name, subsampling=subsamp, celldata=[u_h])
    print("*************************************")
    print("**** Completed simulation",name)
    print("*************************************")
    return u_h
