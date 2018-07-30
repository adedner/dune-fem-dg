from ufl import *
from dune.ufl import Space

def CompressibleEuler(dim, gamma):
    class Model:
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
        def outflowFlux(t,x,u,n):
            return Model.F_c(t,x,u)*n
        boundaryFlux = {}
        for i in range(1,5): boundaryFlux.update( {i: outflowFlux} )
    return Model
def CompressibleEulerDirichlet(dim, gamma):
    class Model(CompressibleEuler(dim,gamma)):
        def outflowValue(t,x,u):
            return u
        boundaryValue = {}
        for i in range(1,5): boundaryValue.update( {i: outflowValue} )
    return Model
def CompressibleEulerSlip(dim, gamma):
    class Model(CompressibleEuler(dim,gamma)):
        def outflowFlux(t,x,u,n):
            _,_, p = CompressibleEuler(dim,gamma).toPrim(u)
            return as_vector([ 0, *(p*n), 0 ])
        boundaryFlux = {}
        for i in range(1,5): boundaryFlux.update( {i: outflowFlux} )
    return Model

def riemanProblem(x,x0,UL,UR):
    return conditional(x<x0,UL,UR)

def constant(dim,gamma):
    return CompressibleEulerDirichlet(dim,gamma) ,\
           as_vector( [0.1,0.,0.,0.1] ),\
           [-1, 0], [1, 0.1], [50, 5],\
           "constant"
def sod(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerNeuman(dim,gamma) ,\
           riemanProblem( x[0], 0.,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-1, 0], [1, 0.1], [50, 5],\
           "sod"
def radialSod1(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerDirichlet(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-1, -1], [1, 1], [50, 50],\
           "radialSod1"
def radialSod2(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerNeuman(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1]),
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1])),\
           [-1, -1], [1, 1], [50, 50],\
           "radialSod2"
def radialSod3(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return CompressibleEulerSlip(dim,gamma) ,\
           riemanProblem( sqrt(dot(x,x)), 0.3,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1])),\
           [-1, -1], [1, 1], [50, 50],\
           "radialSod3"
