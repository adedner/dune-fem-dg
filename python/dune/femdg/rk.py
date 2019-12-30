import math, time, sys
from dune.grid import structuredGrid, cartesianDomain, OutputType
import dune.create as create
from dune.fem.space import dgonb
from dune.fem.function import integrate
from dune.ufl import Constant
from ufl import dot, SpatialCoordinate
from dune.femdg import femDGOperator, rungeKuttaSolver

class FemDGStepper:
    def __init__(self,op,parameters):
        self.rkScheme = rungeKuttaSolver( op, parameters=parameters )
    def __call__(self,u,dt=None):
        if dt is not None:
            self.rkScheme.setTimeStepSize(dt)
        self.rkScheme.solve(u)
        return self.rkScheme.deltaT()
def femdgStepper(parameters):
    def _femdgStepper(op,cfl=None):
        if cfl is not None:
            parameters["fem.timeprovider.factor"] = cfl
        elif not "fem.timeprovider.factor" in parameters:
            parameters["fem.timeprovider.factor"] = 0.45
        return FemDGStepper(op,parameters)
    return _femdgStepper

# http://www.sspsite.org
# https://arxiv.org/pdf/1605.02429.pdf
# https://openaccess.leidenuniv.nl/bitstream/handle/1887/3295/02.pdf?sequence=7
# https://epubs.siam.org/doi/10.1137/07070485X
class Heun:
    def __init__(self, op, cfl=None):
        self.stages = 2
        self.op = op
        self.c = [0,1]
        self.b = [0.5,0.5]
        self.A = [[0,0],[1,0]]
        self.k = self.stages*[None]
        self.cfl = 0.45 if cfl is None else cfl
        for i in range(self.stages):
            self.k[i] = op.space.interpolate(op.space.dimRange*[0],name="k")
        self.tmp = op.space.interpolate(op.space.dimRange*[0],name="tmp")
    def __call__(self,u,dt=None):
        self.op(u,self.k[0])
        if dt is None:
            dt = self.op.timeStepEstimate*self.cfl
        for i in range(1,self.stages):
            self.tmp.assign(u)
            for j in range(i):
                self.tmp.axpy(dt*self.A[i][j],self.k[j])
            self.op.stepTime(self.c[i],dt)
            self.op(self.tmp,self.k[i])
        for i in range(self.stages):
            u.axpy(dt*self.b[i],self.k[i])
        return dt
class ExplSSP3:
    def __init__(self,stages,op,cfl=None):
        self.op     = op
        self.n      = math.ceil(math.sqrt(stages))
        self.stages = self.n*self.n
        self.r      = self.stages-self.n
        self.q2     = op.space.interpolate(op.space.dimRange*[0],name="q2")
        self.tmp    = self.q2.copy()
        self.cfl    = 0.45 * stages*(1-1/self.n) \
                      if cfl is None else cfl
        for i in range(1,self.stages+1):
            print(i,self.c(i))
    def c(self,i):
        return (i-1)/(self.n*self.n-self.n) \
               if i<=(self.n+2)*(self.n-1)/2 \
               else (i-self.n)/(self.n*self.n-self.n)
    def __call__(self,u,dt=None):
        if dt is None:
            self.op(u, self.tmp)
            dt = self.op.timeStepEstimate*self.cfl
        fac = dt/self.r
        i = 1
        while i <= (self.n-1)*(self.n-2)/2:
            self.op.stepTime(self.c(i),dt)
            self.op(u,self.tmp)
            u.axpy(fac, self.tmp)
            i += 1
        self.q2.assign(u)
        while i <= self.n*(self.n+1)/2:
            self.op.stepTime(self.c(i),dt)
            self.op(u,self.tmp)
            u.axpy(fac, self.tmp)
            i += 1
        u.as_numpy[:] *= (self.n-1)/(2*self.n-1)
        u.axpy(self.n/(2*self.n-1), self.q2)
        while i <= self.stages:
            self.op.stepTime(self.c(i),dt)
            self.op(u,self.tmp)
            u.axpy(fac, self.tmp)
            i += 1
        return dt
def explSSP3(stages):
    return lambda op,cfl=None: ExplSSP3(stages,op,cfl)
class ExplSSP4_10:
    def __init__(self, op,cfl=None):
        self.op     = op
        self.stages = 10
        self.q2     = op.space.interpolate(op.space.dimRange*[0],name="q2")
        self.tmp    = self.q2.copy()
        self.cfl    = 0.45 * self.stages*0.6 if cfl is None else cfl
        for i in range(1,self.stages+1):
            print(i,self.c(i))
    def c(self,i):
        return (i-1)/6 if i<=5 else (i-4)/6
    def __call__(self,u,dt=None):
        if dt is None:
            self.op(u, self.tmp)
            dt = self.op.timeStepEstimate*self.cfl # how is this reset for the next time step
        i = 1
        self.q2.assign(u)
        while i <= 5:
            self.op.stepTime(self.c(i),dt)
            self.op(u,self.tmp)
            u.axpy(dt/6, self.tmp)
            i += 1

        self.q2.as_numpy[:] *= 1/25
        self.q2.axpy(9/25,u)
        u.as_numpy[:] *= -5
        u.axpy(15,self.q2)

        while i <= 9:
            self.op.stepTime(self.c(i),dt)
            self.op(u,self.tmp)
            u.axpy(dt/6, self.tmp)
            i += 1

        self.op.stepTime(self.c(i),dt)
        self.op(u,self.tmp)
        u.as_numpy[:] *= 3/5
        u.axpy(1,self.q2)
        u.axpy(dt/10, self.tmp)
        return dt

def run(Model, Stepper,
        initial=None, domain=None, endTime=None, name=None, exact=None, # deprecated - now Model attributes
        polOrder=1, limiter="default", startLevel=0,
        primitive=None, saveStep=None, subsamp=0,
        dt=None,cfl=None,grid="yasp", space="dgonb", threading=True,
        parameters={}):
    if initial is not None:
        print("deprecated usage of run: add initial etc as class atttributes to the Model")
    else:
        initial = Model.initial
        domain = Model.domain
        endTime = Model.endTime
        name = Model.name
        exact = Model.exact
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
    grid.hierarchicalGrid.globalRefine(startLevel)
    dimR     = Model.dimRange
    t        = 0
    tcount   = 0
    saveTime = saveStep

    # create discrete function space
    space = create.space( space, grid, order=polOrder, dimRange=dimR)
    operator = femDGOperator(Model, space, limiter=limiter, threading=True, parameters=parameters )
    stepper  = Stepper(operator, cfl)
    # create and initialize solution
    u_h = space.interpolate(initial, name='u_h')
    operator.applyLimiter( u_h )

    # preparation for output
    if saveStep is not None:
        x = SpatialCoordinate(space.cell())
        tc = Constant(0.0,"time")
        try:
            velo = [create.function("ufl",space.grid, ufl=Model.velocity(tc,x,u_h), order=2, name="velocity")]
        except AttributeError:
            velo = None
        if name is not None:
            vtk = grid.sequencedVTK(name, subsampling=subsamp,
                   celldata=[u_h],
                   pointdata=primitive(Model,u_h) if primitive else [u_h],
                   cellvector=velo
                )
        else:
            vtk = lambda : None
        try:
            velo[0].setConstant("time",[t])
        except:
            pass
    vtk()

    # measure CPU time
    start = time.time()

    fixedDt = dt is not None

    while operator.time < endTime:
        assert not math.isnan( u_h.scalarProductDofs( u_h ) )
        dt = stepper(u_h)
        # check that solution is meaningful
        if math.isnan( u_h.scalarProductDofs( u_h ) ):
            vtk()
            print('ERROR: dofs invalid t =', t,flush=True)
            print('[',tcount,']','dt = ', dt, 'time = ',t, flush=True )
            sys.exit(1)
        # increment time and time step counter
        operator.addToTime(dt)
        tcount += 1
        print('[',tcount,']','dt = ', dt, 'time = ',operator.time,
                'dtEst = ',operator.timeStepEstimate,
                'elements = ',grid.size(0), flush=True )
        if operator.time > saveTime:
            try:
                velo[0].setConstant("time",[operator.time])
                vtk()
            except:
                pass
            saveTime += saveStep
    runTime = time.time()-start
    print("time loop:",runTime,flush=True)
    print("number of time steps ", tcount,flush=True)
    try:
        velo[0].setConstant("time",[operator.time])
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
        u = exact(tc)
        tc.value = operator.time
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
    return u_h, [error, runTime, tcount]
