from ufl import *
from dune.ufl import Space

def CompressibleEuler(dim, gamma):
    class Model:
        def velocity(U):
            return as_vector( [U[i]/U[0] for i in range(1,dim+1)] )
        def rhoeps(U):
            v = Model.velocity(U)
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
            return U[0], Model.velocity(U), Model.pressure(U)
        def F_c(U):
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
        def maxLambda(U,n):
            rho, v, p = Model.toPrim(U)
            return abs(dot(v,n)) + sqrt(gamma*p/rho)
        def physical(U):
            return conditional( (U[0]>1e-8), conditional( Model.rhoeps(U) > 1e-8 , 1, 0 ), 0 )
        def jump(U,V):
            pL = Model.pressure(U)
            pR = Model.pressure(V)
            return (pL - pR)/(0.5*(pL + pR))
    return Model

def CompressibleEulerNeuman(dim, gamma):
    class Model(CompressibleEuler(dim,gamma)):
        def outflowFlux(x,u,n):
            return Model.F_c(u)*n
        boundaryFlux = {1: outflowFlux}
    return Model
def CompressibleEulerDirichlet(dim, gamma):
    class Model(CompressibleEuler(dim,gamma)):
        def outflowValue(x,u):
            return u
        boundaryValue = {1: outflowValue}
    return Model
def CompressibleEulerSlip(dim, gamma):
    class Model(CompressibleEuler(dim,gamma)):
        def outflowFlux(x,u,n):
            _,_, p = CompressibleEuler(dim,gamma).toPrim(u)
            return as_vector([ 0, *(p*n), 0 ])
        boundaryFlux = {1: outflowFlux}
    return Model

def riemanProblem(x,x0,UL,UR):
    return conditional(x<x0,UL,UR)

def constant(dim,gamma):
    return CompressibleEulerDirichlet(dim,gamma) ,\
           as_vector( [0.1,0.,0.,0.1] ),\
           [-1, 0], [1, 0.1], [50, 5]
def sod(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerNeuman(dim,gamma) ,\
           riemanProblem( x[0], 0.,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-1, 0], [1, 0.1], [50, 5]
def radialSod1(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerDirichlet(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-1, -1], [1, 1], [50, 50]
def radialSod2(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerNeuman(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1]),
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1])),\
           [-1, -1], [1, 1], [50, 50]
def radialSod3(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerSlip(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-1, -1], [1, 1], [50, 50]
