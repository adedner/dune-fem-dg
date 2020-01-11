import time, math, sys
from ufl import *
from dune.ufl import DirichletBC
from dune.grid import structuredGrid, cartesianDomain
from dune.alugrid import aluCubeGrid as Grid
from dune.fem.space import lagrange, dgonb, raviartThomas
from dune.fem.scheme import galerkin
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper

eps = 0.01                # diffusion rate
K = 10                    # reaction rate
Q = 1                     # source strength
P1 = as_vector([0.1,0.1]) # midpoint of first source
P2 = as_vector([0.9,0.9]) # midpoint of second source
RF = 0.075                # radius of sources

def problem():
    gridView = Grid(cartesianDomain([0,0],[1,1],[50,50]))

    def computeVelocityCurl():
        # TODO without dimRange the 'vectorization' is not carried out correctly
        streamSpace = lagrange(gridView, order=1, dimRange=1)
        Psi  = streanSpace.interpolate(0,name="streamFunction")
        u    = TrialFunction(streamSpace)
        phi  = TestFunction(streamSpace)
        x    = SpatialCoordinate(streamSpace)
        form = ( inner(grad(u),grad(phi)) -
                 0.1*2*(2*pi)**2*sin(2*pi*x[0])*sin(2*pi*x[1]) * phi[0] ) * dx
        dbc  = DirichletBC(streamSpace,[0])
        streamScheme = galerkin([form == 0, dbc], solver="cg")
        streamScheme.solve(target=Psi)
        velocity = as_vector([-Psi[0].dx(1),Psi[0].dx(0)]) # note: not extra projection
        return velocity

    def computeVelocityGrad():
        # could also use dg here
        pressureSpace = lagrange(gridView, order=1, dimRange=1)
        pressure      = pressureSpace.interpolate(0,name="pressure")
        u    = TrialFunction(pressureSpace)
        phi  = TestFunction(pressureSpace)
        x    = SpatialCoordinate(pressureSpace)
        n    = FacetNormal(pressureSpace)
        form = inner(grad(u),grad(phi)) * dx
        # dbc  = DirichletBC(pressureSpace,[10*(x[0]-0.5)**3*(x[1]-1/2)**3])
        dbc  = DirichletBC(pressureSpace,[ sin(2*pi*(x[0]-0.5)*(x[1]-0.5)) ])
        pressureScheme = galerkin([form == 0, dbc], solver="cg")
        pressureScheme.solve(target=pressure)
        pressure.plot()
        velocity = grad(pressure[0])
        velocitySpace = raviartThomas(gridView,order=1,dimRange=1)
        velo = velocitySpace.project(velocity,name="velocity")
        velo.plot(vectors=[0,1])
        return velocitySpace.project(velocity,name="velocity")

    class Model:
        # transportVelocity = computeVelocityCurl()
        transportVelocity = computeVelocityGrad()
        dimRange = 3
        def S_ns(t,x,U,DU): # or S_ns for a non stiff source
            # sources
            f1 = conditional(dot(x-P1,x-P1) < RF**2, Q, 0)
            f2 = conditional(dot(x-P2,x-P2) < RF**2, Q, 0)
            f  = conditional(t<5, as_vector([f1,f2,0]), as_vector([0,0,0]))
            # reaction rates
            r = K*as_vector([U[0]*U[1], U[0]*U[1], -2*U[0]*U[1]]) #  + U[2]])
            return f - r
        def F_c(t,x,U):
            return as_matrix([ [*(Model.velocity(t,x,u)*u)] for u in U ])
        def maxLambda(t,x,U,n):
            return abs(dot(Model.velocity(t,x,U),n))
        def velocity(t,x,U):
            return Model.transportVelocity
        def F_v(t,x,U,DU):
            return eps*DU
        def maxDiffusion(t,x,U):
           return eps
        def physical(U):
            return conditional(U[0]>=0,1,0)*conditional(U[1]>=0,1,0)*conditional(U[2]>=0,1,0)
        def dirichletValue(t,x,u):
            return as_vector(Model.dimRange*[0])
        boundary = {(1,2,3,4): dirichletValue}

    Model.initial=as_vector([0,0,0])
    Model.domain=gridView
    Model.endTime=10
    Model.name="chemical"
    return Model

Model = problem()
Stepper = femdgStepper(order=3)

# create discrete function space
space = dgonb( Model.domain, order=3, dimRange=Model.dimRange)
operator = femDGOperator(Model, space, limiter="scaling", threading=True)
stepper  = Stepper(operator)
# create and initialize solution
u_h = space.interpolate(Model.initial, name='u_h')
operator.applyLimiter( u_h )

vtk = Model.domain.sequencedVTK("chemicalreaction", subsampling=1, pointdata=[u_h])
vtk() # output initial solution

t        = 0
tcount   = 0
saveStep = Model.endTime/100
saveTime = saveStep
while t < Model.endTime:
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
    if 1: # t > saveTime:
        print('[',tcount,']','dt = ', dt, 'time = ',t,
                'dtEst = ',operator.timeStepEstimate,
                'elements = ',Model.domain.size(0), flush=True )
        vtk()
        saveTime += saveStep

vtk() # output final solution
