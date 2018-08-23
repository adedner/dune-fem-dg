from ufl import *
from dune.ufl import Space

def ShallowWater(topography):
    dim = 2
    space = Space(2,3)
    x = SpatialCoordinate(space.cell())
    g = 9.81
    class Model:
        dimension = dim+1
        def velo(U):
            return as_vector( [U[i]/U[0] for i in range(1,dim+1)] )
        def toCons(U,x=x):
            v = as_vector( [U[i] for i in range(1,dim+1)] )
            return as_vector( [U[0]-topography(x), *(U[0]*v)] )
        def toPrim(U,x=x):
            return U[0]+topography(x), Model.velo(U)
        def F_c(t,x,U):
            assert dim==2
            h = U[0]
            v = Model.velo(U)
            p = 0.5*g*h*h
            return as_matrix([
                    [h*v[0], h*v[1]],
                    [h*v[0]*v[0] + p, h*v[0]*v[1]],
                    [h*v[0]*v[1], h*v[1]*v[1] + p] ] )
        def maxLambda(t,x,U,n):
            h = U[0]
            v = Model.velo(U)
            return abs(dot(v,n)) + sqrt(g*h)
        def velocity(t,x,U):
            return Model.velo(U)
        def physical(U):
            return conditional( (U[0]>1e-8), 1, 0 )
        def jump(U,V):
            hL = U[0]
            hR = V[0]
            return (hL - hR)/(0.5*(hL + hR))
        def S_ns(t,x,U,DU): # or S_s for a stiff source
            return as_vector( [0, *(-U[0]*g*grad(topography(x))) ])
        # boundary = {range(1,5): lambda t,x,u,n: Model.F_c(t,x,u)*n}
        boundary = {range(1,5): lambda t,x,u: u}
    return Model

def leVeque():
    cutoff = lambda x: conditional(2./5.<x[0],1,0)*conditional(x[0]<3./5.,1,0)
    topography = lambda x: 1./4.*(cos(10*pi*(x[0]-0.5))+1)*cutoff(x)
    space = Space(2,3)
    x = SpatialCoordinate(space.cell())
    initial = conditional(x[0]<0.1,1,conditional(x[0]<0.2,1.2,1))
    return ShallowWater(topography),\
           as_vector( [initial,0,0] ),\
           [0, 0], [1, 0.25], [64, 16], 0.1,\
           "leVeque"
