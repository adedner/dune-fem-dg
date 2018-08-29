from ufl import *
from dune.ufl import Space

def CompressibleEuler(dim, gamma):
    class Model:
        dimension = dim+2
        def velo(U):
            return as_vector( [U[i]/U[0] for i in range(1,dim+1)] )
        def rhoeps(U):
            v = Model.velo(U)
            kin = dot(v,v) * 0.5*U[0]
            rE = U[dim+1]
            return rE - kin
        def pressure(U):
            return (gamma-1)*Model.rhoeps(U)
        def toCons(U):
            v = as_vector( [U[i] for i in range(1,dim+1)] )
            rhoEps = U[dim+1]/(gamma-1.)
            kin = dot(v,v) * 0.5*U[0]
            return as_vector( [U[0], *(U[0]*v), rhoEps+kin] )
        def toPrim(U):
            return U[0], Model.velo(U), Model.pressure(U)

        # interface methods
        def F_c(t,x,U):
            assert dim==2
            rho, v, p = Model.toPrim(U)
            rE = U[dim+1]
            # TODO make indpendent of dim
            res = as_matrix([
                    [rho*v[0], rho*v[1]],
                    [rho*v[0]*v[0] + p, rho*v[0]*v[1]],
                    [rho*v[0]*v[1], rho*v[1]*v[1] + p],
                    [(rE+p)*v[0], (rE+p)*v[1]] ] )
            return res
        def maxLambda(t,x,U,n):
            rho, v, p = Model.toPrim(U)
            return abs(dot(v,n)) + sqrt(gamma*p/rho)
        def velocity(t,x,U):
            return Model.velo(U)
        def physical(U):
            return conditional( (U[0]>1e-8), conditional( Model.rhoeps(U) > 1e-8 , 1, 0 ), 0 )
        def jump(U,V):
            pL = Model.pressure(U)
            pR = Model.pressure(V)
            return (pL - pR)/(0.5*(pL + pR))
    return Model

def CompressibleEulerNeuman(dim, gamma):
    class Model(CompressibleEuler(dim,gamma)):
        boundary = {range(1,5): lambda t,x,u,n: Model.F_c(t,x,u)*n}
    return Model
def CompressibleEulerDirichlet(dim, gamma):
    class Model(CompressibleEuler(dim,gamma)):
        boundary = {range(1,5): lambda t,x,u: u}
    return Model
def CompressibleEulerSlip(dim, gamma):
    class Model(CompressibleEuler(dim,gamma)):
        def outflowFlux(t,x,u,n):
            _,_, p = CompressibleEuler(dim,gamma).toPrim(u)
            return as_vector([ 0, *(p*n), 0 ])
        boundary = {range(1,5): outflowFlux}
    return Model

def riemanProblem(x,x0,UL,UR):
    return conditional(x<x0,UL,UR)

# TODO Add exact solution where available (last argument)
def constant(dim,gamma):
    return CompressibleEulerDirichlet(dim,gamma) ,\
           as_vector( [0.1,0.,0.,0.1] ),\
           [-1, 0], [1, 0.1], [50, 5], 0.1,\
           "constant", None
def sod(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerDirichlet(dim,gamma) ,\
           riemanProblem( x[0], 0.5,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [0, 0], [1, 0.25], [64, 16], 0.15,\
           "sod", None
def radialSod1(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerDirichlet(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-0.5, -0.5], [0.5, 0.5], [20, 20], 0.25,\
           "radialSod1", None
def radialSod1Large(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerDirichlet(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-1.5, -1.5], [1.5, 1.5], [60, 60], 0.5,\
           "radialSod1Large", None
def radialSod2(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerNeuman(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1]),
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1])),\
           [-0.5, -0.5], [0.5, 0.5], [20, 20], 0.25,\
           "radialSod2", None
def radialSod3(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerSlip(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-0.5, -0.5], [0.5, 0.5], [20, 20], 0.5,\
           "radialSod3", None

def leVeque(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    initial = conditional(abs(x[0]-0.15)<0.05,1.2,1)
    return CompressibleEulerDirichlet(dim,gamma),\
           as_vector( [initial,0,0,initial] ),\
           [0, 0], [1, 0.25], [64, 16], 0.7,\
           "leVeque1D", None

def vortex(dim,gamma):
    S = 13.5    # strength of vortex
    R = 1.5     # radius of vortex
    M = 0.4     # Machnumber

    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    f = (1 - x[0]*x[0] - x[1]*x[1])/(2*R*R)
    rho = pow(1 - S*S*M*M*(gamma - 1)*exp(2*f)/(8*pi*pi), 1/(gamma - 1))
    u =     S*x[1]*exp(f)/(2*pi*R)
    v = 1 - S*x[0]*exp(f)/(2*pi*R)
    p = rho / (gamma*M*M)
    return CompressibleEulerDirichlet(dim,gamma),\
           as_vector( [rho,u,v,p] ),\
           [-2, 2], [2, 2], [64, 64], 100,\
           "vortex", None
