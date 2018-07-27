from ufl import *
from dune.ufl import Space
class Transport1D:
    name = 'transport'

    def F_c(U):
        return [ [U[0], 0] ]

    def outflowValue(u):
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

class Burgers1D:
    name = 'burgers'

    def F_c(U):
        return [ [U[0]*U[0]/2, 0] ]

    def outflowValue(u):
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
def RiemanProblem(UL,UR):
    x = SpatialCoordinate(space.cell())
    return conditional(x[0]<-0.5,as_vector(UL),as_vector(UR))
def Constant():
    return as_vector( [0.1] )
