from ufl import *
from dune.ufl import Space
class Transport1D:
    name = 'transport'

    def F_c(U):
        return [ [U[0], 0] ]

    def outflowValue(x,u):
        return u
    boundaryValue = {1: outflowValue}

    def maxLambda(U,n):
        return abs(n[0])
    def velocity(U):
        return [1,0]
    def physical(U):
        return 1
    def jump(U,V):
        return U[0]-V[0]

class LinearAdvectionDiffusion1D:
    name = 'linearAD'

    def F_c(U):
        return [ [U[0], 0] ]
    def F_v(U,DU):
        return 0.1*DU

    def dirichletValue(x,u):
        return as_vector( [cos(2*pi*x[0])] )
    boundaryValue = {1: dirichletValue}

    def maxLambda(U,n):
        return abs(n[0])
    def velocity(U):
        return [1,0]
    def physical(U):
        return 1
    def jump(U,V):
        return U[0]-V[0]

class Burgers1D:
    name = 'burgers'

    def F_c(U):
        return [ [U[0]*U[0]/2, 0] ]

    def outflowValue(x,u):
        return u
    boundaryValue = {1: outflowValue}

    def maxLambda(U,n):
        return abs(U[0]*n[0])
    def velocity(U):
        return [U[0],0]
    def physical(U):
        return 1
    def jump(U,V):
        return U[0]-V[0]

space = Space(2,1)
x = SpatialCoordinate(space.cell())
def riemanProblem(x,x0,UL,UR):
    return conditional(x<x0,as_vector(UL),as_vector(UR))

def constantTransport():
    return Transport1D, as_vector( [0.1] )

def cosProblem():
    return LinearAdvectionDiffusion1D, as_vector( [cos(2*pi*x[0])] )

def burgersShock():
    return Burgers1D, riemanProblem(x[0],-0.5,[1],[0])

def burgersVW():
    return Burgers1D, riemanProblem(x[0],0,[-1],[1])

def burgersStationaryShock():
    return Burgers1D, riemanProblem(x[0],0,[1],[-1])
