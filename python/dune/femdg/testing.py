import time, math, sys
from dune.grid import structuredGrid, cartesianDomain, OutputType
import dune.create as create
from dune.fem.function import integrate
from dune.ufl import Constant, expression2GF
from ufl import dot, SpatialCoordinate
from dune.femdg import femDGOperator, rungeKuttaSolver, createLimiter
from dune.femdg.rk import femdgStepper

from collections import namedtuple

Parameters = namedtuple("TestingParameters",
                        ["Model", "initial", "domain", "endTime", "name", "exact"])
Parameters.__new__.__defaults__ = (None,None) # defaults for name, exact


def run(Model, Stepper=None,
        initial=None, domain=None, endTime=None, name=None, exact=None, # deprecated - now Model attributes
        polOrder=1, limiter="default", startLevel=0,
        primitive=None, saveStep=None, subsamp=0,
        dt=None,cfl=None,grid="yasp", space="dgonb", threading=True,
        parameters={}):
    if Stepper is None:
        Stepper = femdgStepper(parameters=parameters)
    if initial is not None:
        print("deprecated usage of run: add initial etc as class atttributes to the Model")
    else:
        initial = Model.initial
        domain = Model.domain
        endTime = Model.endTime
        name = Model.name
        try:
            exact = Model.exact
        except:
            exact = None
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
        # periodic boundaries do not work for YaspGrid because the concept
        # is different from ALUGrid and SPGrid which is the concept used in
        # dune-fem-dg
        elif grid == 'yasp':
           errorstr = "YaspGrid (grid='"+grid+"') does not work with periodic  boundaries!\nChoose grid='alucube' or grid='spbisection'"
           raise Exception(errorstr)

        # create domain and grid
        domain   = cartesianDomain(x0,x1,N,periodic=periodic,overlap=0)
        grid     = create.grid(grid,domain)
        print("Setting periodic boundaries",periodic,flush=True)
    except TypeError: # assume the 'domain' is already a gridview
        grid = domain
    # initial refinement of grid
    if startLevel>0: # note this call clears all discrete function e.g. the velocity in the chemical problem
        grid.hierarchicalGrid.globalRefine(startLevel)
    dimR     = Model.dimRange
    t        = 0
    tcount   = 0
    saveTime = saveStep

    # create discrete function space
    space = create.space( space, grid, order=polOrder, dimRange=dimR)
    operator = femDGOperator(Model, space, limiter=limiter, threading=threading, parameters=parameters )
    stepper  = Stepper(operator, cfl)
    # create and initialize solution
    u_h = space.interpolate(initial, name='u_h')
    operator.applyLimiter( u_h )

    # preparation for output
    vtk = lambda : None
    if saveStep is not None and name is not None:
        x = SpatialCoordinate(space.cell())
        tc = Constant(0.0,"time")
        try:
            velo = create.function("ufl",space.grid, ufl=Model.velocity(tc,x,u_h), order=2, name="velocity")
        except AttributeError:
            velo = None
        vtk = grid.sequencedVTK(name, subsampling=subsamp,
               celldata=[u_h],
               pointdata=primitive(Model,u_h) if primitive else [u_h],
               cellvector=[velo]
            )
        try:
            velo.setConstant("time",[t])
        except:
            pass
    vtk()

    # measure CPU time
    start = time.time()

    fixedDt = dt is not None

    # tracemalloc.start()

    t = 0
    while t < endTime:
        assert not math.isnan( u_h.scalarProductDofs( u_h ) )
        operator.setTime(t)
        dt = stepper(u_h)
        t += dt
        tcount += 1
        # check that solution is meaningful
        if math.isnan( u_h.scalarProductDofs( u_h ) ):
            vtk()
            print('ERROR: dofs invalid t =', t,flush=True)
            print('[',tcount,']','dt = ', dt, 'time = ',t, flush=True )
            sys.exit(1)
        if t > saveTime:
            print('[',tcount,']','dt = ', dt, 'time = ',t,
                    'dtEst = ',operator.timeStepEstimate,
                    'elements = ',grid.size(0), flush=True )
            try:
                velo0.setConstant("time",[t])
            except:
                pass
            vtk()
            saveTime += saveStep
        # snapshot = tracemalloc.take_snapshot()
        # display_top(snapshot)

    runTime = time.time()-start
    print("time loop:",runTime,flush=True)
    print("number of time steps ", tcount,flush=True)
    try:
        velo.setConstant("time",[t])
        vtk()
    except:
        pass

    # output the final result and compute error (if exact is available)
    if exact is not None and name is not None:
        tc = Constant(0, "time")
        # using '0' here because uusing 't'
        # instead means that this value is added to the generated code so
        # that slight changes to 't' require building new local functions -
        # should be fixed in dune-fem
        u = expression2GF(grid,exact(tc),order=5)
        tc.value = t
        grid.writeVTK(name, subsampling=subsamp,
                celldata=[u_h], pointdata={"exact":u})
        error = integrate( grid, dot(u-u_h,u-u_h), order=5 )
    elif name is not None:
        grid.writeVTK(name, subsampling=subsamp, celldata=[u_h])
        error = integrate( grid, dot(u_h,u_h), order=5 )
    error = math.sqrt(error)
    print("*************************************")
    print("**** Completed simulation",name)
    print("**** error:", error)
    print("*************************************",flush=True)
    return u_h, [error, runTime, tcount, operator.counter()]


def oldRun(Model,
        initial=None, domain=None, endTime=None, name=None, exact=None, # deprecated - now Model attributes
        polOrder=1, limiter="default", startLevel=0,
        primitive=None, saveStep=None, subsamp=0,
        dt=None,grid="yasp", space="dgonb", threading=True,
        parameters={}):
    if initial is not None:
        print("deprecated usage of run: add initial etc as class atttributes to the Model")
    else:
        initial = Model.initial
        domain = Model.domain
        endTime = Model.endTime
        name = Model.name
        try:
            exact = Model.exact
        except:
            exact = None
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
        # create domain and grid
        domain   = cartesianDomain(x0,x1,N,periodic=periodic,overlap=0)
        grid     = create.grid(grid,domain)
        print("Setting periodic boundaries",periodic,flush=True)
    except TypeError: # assume the 'domain' is already a gridview
        grid = domain
    # initial refinement of grid
    if startLevel>0: # note this call clears all discrete function e.g. the velocity in the chemical problem
        grid.hierarchicalGrid.globalRefine(startLevel)

    dimR     = Model.dimRange
    t        = 0
    count    = 0
    saveTime = saveStep

    # create discrete function space
    space = create.space( space, grid, order=polOrder, dimRange=dimR)
    # create and initialize solution
    u_h = space.interpolate(initial, name='u_h')
    u_h_n = space.interpolate(initial, name='u_h_n')
    # create solution scheme, i.e. operator and ODE solver

    scalingLimiter=None
    if limiter == "scaling":
        limiter = None
        scalingLimiter = createLimiter( space, limiter='scaling' )

    # create DG operator based on Model
    operator = femDGOperator(Model, space, limiter=limiter, threading=threading, parameters=parameters )

    # create Runge-Kutta solver
    rkScheme = rungeKuttaSolver( operator, parameters=parameters )

    # limit initial data if necessary
    operator.applyLimiter( u_h )

    print("number of elements: ",grid.size(0),flush=True)

    # preparation for output
    if saveStep is not None:
        x = SpatialCoordinate(space.cell())
        tc = Constant(0.0,"time")
        try:
            velo = create.function("ufl",space.grid, ufl=Model.velocity(tc,x,u_h), order=2, name="velocity")
        except AttributeError:
            velo = None
        # BUG: the local function should have a plot command - like the df
        # velo.plot()
        if name is not None:
            vtk = grid.writeVTK(name, subsampling=subsamp, write=False,
                   celldata=[u_h],
                   pointdata=primitive(Model,u_h) if primitive else [u_h],
                   cellvector=[velo]
                )
        else:
            vtk = None
        try:
            velo.setConstant("time",[t])
        except:
            pass
        if vtk is not None:
            vtk.write(name, count) # , outputType=OutputType.appendedraw)

    # measure CPU time
    start = time.time()

    tcount = 0
    # time loop
    # set fixed time step size to ODE solver
    if dt is not None:
        rkScheme.setTimeStepSize(dt)
    print("Start solving")
    while t < endTime:
        # solver time step
        assert not math.isnan( u_h.scalarProductDofs( u_h ) )
        rkScheme.solve(u_h)
        if scalingLimiter is not None:
            u_h_n.assign( u_h )
            scalingLimiter( u_h_n, u_h )

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
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'elements = ',grid.size(0), 'count = ',count, flush=True )
        if saveStep is not None and t > saveTime:
            count += 1
            try:
                velo.setConstant("time",[t])
            except:
                pass
            if vtk is not None:
                vtk.write(name, count, outputType=OutputType.appendedraw)
            saveTime += saveStep

    print("time loop:",time.time()-start,flush=True)
    print("number of time steps ", tcount,flush=True)
    if name is not None and saveStep is not None:
        try:
            velo.setConstant("time",[t])
        except:
            pass
        vtk.write(name, count, outputType=OutputType.appendedraw)

    # output the final result and compute error (if exact is available)
    if exact is not None and name is not None:
        tc = Constant(0, "time")
        # using '0' here because uusing 't'
        # instead means that this value is added to the generated code so
        # that slight changes to 't' require building new local functions -
        # should be fixed in dune-fem
        u = exact(tc)
        tc.value = t
        grid.writeVTK(name, subsampling=subsamp,
                celldata=[u_h], pointdata={"exact":u})
        error = integrate( grid, dot(u_h-u,u_h-u), order=5 )
    elif name is not None:
        grid.writeVTK(name, subsampling=subsamp, celldata=[u_h])
        error = integrate( grid, dot(u_h,u_h), order=5 )
    error = math.sqrt(error)
    print("*************************************")
    print("**** Completed simulation",name)
    print("**** error:", error)
    print("*************************************",flush=True)
    return u_h, (error,operator.counter)
