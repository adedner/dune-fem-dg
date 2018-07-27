from ufl import *
exec(open("euler.ufl").read())

def RiemanProblem(UL,UR):
    x = SpatialCoordinate(space.cell())
    return conditional(x[0]<0,as_vector(UL),as_vector(UR))
def Sod():
    return RiemanProblem( Model.toCons([1,0,0,1]),
                          Model.toCons([0.125,0,0,0.1]) )
def Constant():
    return as_vector( [0.1,0.,0.,0.1] )
