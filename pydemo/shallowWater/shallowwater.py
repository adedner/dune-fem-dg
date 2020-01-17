from ufl import *
from dune.ufl import Space
from dune.femdg.testing import Parameters
# Parameters are: "Model", "initial", "domain", "endTime", "name", "exact"

def ShallowWater(topography,g):
    dim = 2
    space = Space(2,3)
    x = SpatialCoordinate(space.cell())
    class Model:
        def velo(U):
            return as_vector( [U[i]/U[0] for i in range(1,dim+1)] )
        def toCons(U,x=x):
            return as_vector( [U[0]-topography(x)]+[U[i] for i in range(1,dim+1)] )
        def toPrim(U,x=x):
            return U[0]+topography(x), Model.velo(U)

        # interface methods
        dimRange = dim+1
        def F_c(t,x,U):
            assert dim==2
            h, v = U[0], Model.velo(U)
            p = 0.5*g*h*h
            return as_matrix([
                    [h*v[0], h*v[1]],
                    [h*v[0]*v[0] + p, h*v[0]*v[1]],
                    [h*v[0]*v[1], h*v[1]*v[1] + p] ] )
        def maxLambda(t,x,U,n):
            h,v = U[0], Model.velo(U)
            return abs(dot(v,n)) + sqrt(g*h)
        def velocity(t,x,U):
            return Model.velo(U)
        def physical(U):
            return conditional( U[0]>1e-8, 1, 0 )
        def jump(U,V):
            hL, hR = U[0], V[0]
            return (hL - hR)/(0.5*(hL + hR))
        def S_expl(t,x,U,DU): # or S_impl for a stiff source
            return as_vector( [0, *(-U[0]*g*grad(topography(x))) ])
        boundary = {range(1,5): lambda t,x,u,n: Model.F_c(t,x,u)*n}
        # boundary = {range(1,5): lambda t,x,u: u}

    return Model

# example 5.1 and 7.1 from
# https://www.sciencedirect.com/science/article/pii/S0021999198960582
def leVeque(dim):
    space = Space(2,3)
    x = SpatialCoordinate(space.cell())
    if dim == 1:
        topography = lambda x: conditional(abs(x[0]-0.5)<0.1, 1./4.*(cos(10*pi*(x[0]-0.5))+1), 0)
        initial = conditional(abs(x[0]-0.15)<0.05,1.2,1)
        Model = ShallowWater(topography,1)
        return Parameters(Model=Model,
                          initial=Model.toCons(as_vector( [initial,0,0] )),
                          domain=[[0, 0], [1, 0.25], [64, 16]], endTime=0.7,
                          name="leVeque1D")
    else:
        topography = lambda x: 0.8*exp(-5*(x[0]-0.9)**2-50*(x[1]-0.5)**2)
        initial = conditional(abs(x[0]-0.1)<0.05,1.01,1)
        Model = ShallowWater(topography,1)
        return Parameters(Model=Model,
                          initial=Model.toCons(as_vector( [initial,0,0] )),
                          domain=[[0, 0], [2, 1], [80, 40]], endTime=1.8,
                          name="leVeque2D")

leVeque1D = lambda: leVeque(1)
leVeque2D = lambda: leVeque(2)
problems = [leVeque1D,leVeque2D]
