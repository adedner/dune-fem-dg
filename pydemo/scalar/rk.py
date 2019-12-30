import time, math, sys
from dune.grid import structuredGrid
from dune.fem.space import dgonb
from dune.fem.function import integrate
from dune.ufl import Constant
from ufl import dot
from dune.femdg import femDGOperator

from dune.fem import parameter

# from scalar import shockTransport as problem
# from scalar import sinProblem as problem
# from scalar import sinTransportProblem as problem
# from scalar import sinAdvDiffProblem as problem
from scalar import pulse as problem
# from scalar import diffusivePulse as problem

parameter.append({"fem.verboserank": 0})
parameters = {"fem.ode.odesolver": "EX",   # EX, IM, IMEX
              "fem.ode.order": 3,
              "fem.ode.verbose": "cfl",      # none, cfl, full
              "fem.ode.cflMax": 100,
              "fem.ode.cflincrease": 1.25,
              "fem.ode.miniterations": 35,
              "fem.ode.maxiterations": 100,
              "fem.timeprovider.factor": 0.45,
              "dgadvectionflux.method": "LLF",
              "fem.solver.gmres.restart": 50,
              "dgdiffusionflux.method": "CDG2",      # CDG2, CDG, BR2, IP, NIPG, BO
              "dgdiffusionflux.theoryparameters": 1, # scaling with theory parameters
              "dgdiffusionflux.penalty": 0,
              "dgdiffusionflux.liftfactor": 1}

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
    def __init__(self, stages,op,cfl=None):
        self.op     = op
        self.n      = math.ceil(math.sqrt(stages))
        self.stages = self.n*self.n
        self.r      = self.stages-self.n
        self.q2     = op.space.interpolate(space.dimRange*[0],name="q2")
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
class ExplSSP4_10:
    def __init__(self, op,cfl=None):
        self.op     = op
        self.stages = 10
        self.q2     = op.space.interpolate(space.dimRange*[0],name="q2")
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

def run(Model, space,
        dt=None, cfl=None, limiter="default",
        saveStep=None, subsamp=0):
    endTime  = Model.endTime
    name     = Model.name
    tcount   = 0
    saveTime = saveStep

    stages = 4
    operator = femDGOperator(Model, space, limiter=limiter, threading=True, parameters=parameters )
    # stepper  = Heun(operator, cfl)
    # stepper  = ExplSSP3(stages,operator,cfl)
    stepper  = ExplSSP4_10(operator,cfl)

    u_h = space.interpolate(Model.initial, name='u_h')
    operator.applyLimiter( u_h )
    vtk = grid.sequencedVTK(name, pointdata=[u_h], subsampling=subsamp)
    vtk()
    fixedDt = dt is not None

    while operator.time < endTime:
        dt = stepper(u_h)
        # check that solution is meaningful
        assert not math.isnan( u_h.scalarProductDofs( u_h ) )
        # increment time and time step counter
        operator.addToTime(dt)
        tcount += 1
        print('[',tcount,']','dt = ', dt, 'time = ',operator.time,
                'dtEst = ',operator.timeStepEstimate,
                'elements = ',grid.size(0), flush=True )
        if saveStep is not None and operator.time > saveTime:
            vtk()
            saveTime += saveStep
    print("number of time steps ", tcount,flush=True)
    vtk()
    return u_h, operator.time

Model = problem()
grid  = structuredGrid(*Model.domain)
grid.hierarchicalGrid.globalRefine(2)
space = dgonb(grid, order=3, dimRange=Model.dimRange)

fixedDt = None
u_h, endTime = run(Model, space, limiter=None, dt=fixedDt, saveStep=0.01)

exact = Model.exact
if exact is not None:
    tc = Constant(0, "time")
    u  = exact(tc)
    tc.value = endTime # there is an issue here with the generated code in the next step involving the actual endTime
    error = integrate( grid, dot(u_h-u,u_h-u), order=5 )
    error = math.sqrt(error)
    print("error:", error)
