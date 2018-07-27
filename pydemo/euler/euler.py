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
            return [U[0], *(U[0]*v), rhoEps+kin]
        def toPrim(U):
            return [U[0], Model.velocity(U), Model.pressure(U)]
        def F_c(U):
            assert dim==2
            rho, v, p = Model.toPrim(U)
            rE = U[dim+1]
            res = [ [rho*v[0], rho*v[1]],
                    [rho*v[0]*v[0] + p, rho*v[0]*v[1]],
                    [rho*v[0]*v[1], rho*v[1]*v[1] + p],
                    [(rE+p)*v[0], (rE+p)*v[1]] ]
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

        def outflowFlux(u,n):
            return Model.F_c(u)*n
        def outflowValue(u):
            return u

        # boundaryFlux = {1: outflowFlux}
        boundaryValue = {1: outflowValue}
    return Model

def RiemanProblem(x,x0,UL,UR):
    return conditional(x[0]<x0,as_vector(UL),as_vector(UR))
def Sod(dim,gamma):
    space = Space(dim,dim+2)
    x = SpatialCoordinate(space.cell())
    return RiemanProblem( x, 0.,
                          CompressibleEuler(dim,gamma).toCons([1,0,0,1]),
                          CompressibleEuler(dim,gamma).toCons([0.125,0,0,0.1]) )
def Constant():
    return as_vector( [0.1,0.,0.,0.1] )
