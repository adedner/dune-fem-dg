from ufl import *
from dune.ufl import Space
class Transport1D:
    name = 'transport'

    def F_c(t,x,U):
        return as_matrix( [ [U[0], 0] ] )

    def outflowValue(t,x,u):
        return u
    boundaryValue = {}
    for i in range(1,5): boundaryValue.update( {i: outflowValue} )

    def maxLambda(t,x,U,n):
        return abs(n[0])
    def velocity(t,x,U):
        return as_vector( [1,0] )
    def physical(U):
        return 1
    def jump(U,V):
        return U-V

class LinearAdvectionDiffusion1D:
    name = 'linearAD'

    # def F_c(t,x,U):
    #     return as_matrix( [ [U[0], 0] ] )
    def F_v(t,x,U,DU):
        return 0.5*DU
    def maxLambda(t,x,U,n):
        return abs(n[0])
    def velocity(t,x,U):
        return as_vector( [0,0] )
    def physical(U):
        return 1
    def jump(U,V):
        return U-V
def LinearAdvectionDiffusion1DDirichlet(bnd):
    class Model(LinearAdvectionDiffusion1D):
        def dirichletValue(t,x,u):
            return bnd
        boundaryValue = {}
        for i in range(1,5): boundaryValue.update( {i: dirichletValue} )
    return Model
class LinearAdvectionDiffusion1DNeuman(LinearAdvectionDiffusion1D):
    def zeroFlux(t,x,u,n):
        return 0
    boundaryFlux = {}
    for i in range(1,5): boundaryFlux.update( {i: zeroFlux} )

class Burgers1D:
    name = 'burgers'

    def F_c(t,x,U):
        return as_matrix( [ [U[0]*U[0]/2, 0] ] )

    def outflowValue(t,x,u):
        return u
    boundaryValue = {}
    for i in range(1,5): boundaryValue.update( {i: outflowValue} )

    def maxLambda(t,x,U,n):
        return abs(U[0]*n[0])
    def velocity(t,x,U):
        return as_vector( [U[0],0] )
    def physical(U):
        return 1
    def jump(U,V):
        return U-V

space = Space(2,1)
x = SpatialCoordinate(space.cell())
def riemanProblem(x,x0,UL,UR):
    return conditional(x<x0,as_vector(UL),as_vector(UR))

def constantTransport():
    return Transport1D, as_vector( [0.1] )
def shockTransport():
    return Transport1D, riemanProblem(x[0],-0.5, [1], [0])

def sinProblem():
    u0 = as_vector( [sin(2*pi*x[0])] )
    return LinearAdvectionDiffusion1DDirichlet(u0), u0

def cosProblem():
    return LinearAdvectionDiffusion1DNeuman, as_vector( [cos(2*pi*x[0])] )

def burgersShock():
    return Burgers1D, riemanProblem(x[0],-0.5,[1],[0])

def burgersVW():
    return Burgers1D, riemanProblem(x[0],0,[-1],[1])

def burgersStationaryShock():
    return Burgers1D, riemanProblem(x[0],0,[1],[-1])
