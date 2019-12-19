from ufl import *
from dune.ufl import Space, DirichletBC
from dune.grid import cartesianDomain
from dune.alugrid import aluCubeGrid as Grid
import dune.fem
from dune.fem import parameter
from dune.femdg.testing import run, Parameters


def problem():
    gridView = Grid(cartesianDomain([0,0],[1,1],[50,50]))
    def velocity():
        # TODO without dimRange the 'vectorization' is not carried out correctly
        velocitySpace = dune.fem.space.lagrange(gridView, order=1, dimRange=1)
        Psi  = velocitySpace.interpolate(0,name="streamFunction")
        u    = TrialFunction(velocitySpace)
        phi  = TestFunction(velocitySpace)
        x    = SpatialCoordinate(velocitySpace)
        form = ( inner(grad(u),grad(phi)) -
                 0.1*2*(2*pi)**2*sin(2*pi*x[0])*sin(2*pi*x[1]) * phi[0] ) * dx
        dbc  = DirichletBC(velocitySpace,[0])
        velocityScheme = dune.fem.scheme.galerkin([form == 0, dbc], solver="cg")
        velocityScheme.solve(target=Psi)
        return as_vector([-Psi[0].dx(1),Psi[0].dx(0)])

    v = velocity()

    eps = 0.01                # diffusion rate
    K = 10                    # reaction rate
    Q = 0.1                   # source strength
    P1 = as_vector([0.1,0.1]) # midpoint of first source
    P2 = as_vector([0.9,0.9]) # midpoint of second source
    RF = 0.075                # radius of sources

    class Model:
        dimRange = 3
        def S_ns(t,x,U,DU): # or S_ns for a non stiff source
            # sources
            f1 = conditional(dot(x-P1,x-P1) < RF**2, Q, 0)
            f2 = conditional(dot(x-P2,x-P2) < RF**2, Q, 0)
            f  = as_vector([f1,f2,0])
            # reaction rates
            r = K*as_vector([U[0]*U[1], U[0]*U[1], -2*U[0]*U[1] + 10*U[2]])
            return f - r
        def F_c(t,x,U):
            return as_matrix([ [*(v*u)] for u in U ])
        def maxLambda(t,x,U,n):
            return abs(dot(v,n))
        def velocity(t,x,U):
            return v
        def F_v(t,x,U,DU):
            return eps*DU
        def maxDiffusion(t,x,U):
           return eps
        def physical(U):
            return 1
        def jump(U,V):
            return abs(U-V)
        def dirichletValue(t,x,u):
            return as_vector(Model.dimRange*[0])
        boundary = {(1,2,3,4): dirichletValue}

    return Parameters(Model=Model, initial=as_vector([0,0,0]),
                      domain=gridView, endTime=5, name="chemical")

parameter.append({"fem.verboserank": 0})

parameters = {"fem.ode.odesolver": "IMEX",    # EX, IM, IMEX
              "fem.ode.order": 1,
              "fem.ode.verbose": "none",      # none, cfl, full
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

run(*problem(), dt=0.01,
    startLevel=0, polOrder=1, limiter=None,
    primitive=None, saveStep=0.01, subsamp=0,
    parameters=parameters)
