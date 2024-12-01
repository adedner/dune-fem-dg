import dune.femdg

import numpy as np
from matplotlib import pyplot
from ufl import *
from dune.ufl import DirichletBC, Space
from dune.grid import structuredGrid
from dune.fem.function import gridFunction
from dune.fem.plotting import plotPointData as plot
from dune.fem.space import lagrange, dgonb
from dune.fem.scheme import galerkin
from dune.femdg import femDGUFL
from dune.femdg.rk import femdgStepper

from dune.fem import threading
threading.use = max(4,threading.max) # use at most 4 threads

gridView = structuredGrid([0,0],[1.2*np.pi,2.2*np.pi],[30,60])
coord = SpatialCoordinate(Space(2))

def test(model,Psi, show=False, tol=0):
    if Psi is None:
        space = lagrange(gridView, order=2, dimRange=1)
    else:
        space = Psi.space
    u = space.interpolate(0,name="streamFunction")
    # dune.generator.setNoDependencyCheck()
    scheme = galerkin( femDGUFL(model(Psi), space) )
    # dune.generator.setDependencyCheck()
    scheme.solve(target=u)
    if Psi is not None:
        if show:
            fig,(ax1,ax2,ax3) = pyplot.subplots(1,3,figsize=(20,5))
            clim = None # [-1,1]
            plot(Psi, gridView=gridView, gridLines=None,
                      clim=clim,
                      figure=(fig,ax1))
            plot(u, gridView=gridView, gridLines=None,
                    clim=clim,
                    figure=(fig,ax2))
            plot(Psi-u, gridView=gridView, gridLines=None, figure=(fig,ax3))
            pyplot.show()

        if tol > 0:
            if np.max(Psi.as_numpy-u.as_numpy) >= tol:
                print(tol,np.max(Psi.as_numpy-u.as_numpy))
            assert np.max(Psi.as_numpy-u.as_numpy) < tol
        if tol == 0:
            if not np.all(np.isclose( Psi.as_numpy-u.as_numpy, 0, atol=1e-6) ):
                print(tol,np.max(Psi.as_numpy-u.as_numpy))
            assert np.all(np.isclose( Psi.as_numpy-u.as_numpy, 0, atol=1e-6) )
    return u

# L[u] = L_e[u] + L_i[u] = (-div(F_c) - S_e) + (-div(F_v) - S_i)
# If       F_v = grad(u) then L[u] = -laplace(u) - S_i so S_i = -laplace(u)
# and with F_c = bu then L[u] = -div(bu) - S_e - laplace(u) - S_i
# Weak form: <L[u],v> = (F_c+F_v).grad(v) - (S_e+S_i)v - (F_c+F_v).n v
class BaseModel:
    exact = lambda x: as_vector([ sin(x[0])*cos(x[1]) ])
    b = lambda x: as_vector([ sin(x[0])*sin(x[1]), cos(x[0])*cos(x[1]) ])
    # Dexact   = (cos(x)cos(y), -sin(x)sin(y))
    # lapexact = -sin(x)cos(y) - sin(x)cos(y) = -2exact
    # b.Dexact = sin(x)sin(y) cos(x)cos(y) - cos(x)cos(y) sin(x)sin(y) #
    #          = cos(x)sin(y) exact - cos(x)sin(y) exact
    # top:    exact, Dexact.n = sin(x), -sin(x)sin(y) = 0  (y=2pi)
    # bottom: exact, Dexact.n = sin(x), sin(x)sin(y)  = 0  (y=0)

    def S_i(t,x,U,DU):
        return 2 * BaseModel.exact(x)
    def F_v(t,x,U,DU):
        return DU
    def S_e(t,x,U,DU):
        return div(BaseModel.b(x))*U                # sign?
    def F_c(t,x,U):
        return as_vector([ BaseModel.b(x)*U[0] ])   # sign?

def dirichletTest(Psi):
    if Psi is None: # use the exact solution as Dirichlet boundary conditions
        class Model(BaseModel):
            boundary = {range(1,5): BaseModel.exact(coord)}
    else: # use discrete solution as Dirichlet boundary conditions
        class Model(BaseModel):
            boundary = {range(1,5): Psi}
    return Model

def neumannTest(Psi,exact):
    class Model(BaseModel):
        def bndFlux_c(t,x,U,n):
            if exact:
                return dot(Model.F_c(t,x,BaseModel.exact(x)),n)
                # - dot(Model.b(x),n) * BaseModel.exact(x)  # sign?
            else:
                return dot(Model.F_c(t,x,Psi),n)
                # - dot(Model.b(x),n) * Psi                 # sign?
        def bndFlux_v(t,x,U,DU,n):
            if exact:
                return dot(Model.F_v(t,x,BaseModel.exact(x),grad(BaseModel.exact(x))),n)
                # return as_vector([ dot(grad(BaseModel.exact(x)[0]),n) ])
            else:
                return dot(Model.F_v(t,x,Psi,grad(Psi)),n)
                # return as_vector([ dot(grad(Psi[0]),n) ])
        if exact:
            boundary = {(1,2): BaseModel.exact(coord), (3,4): (bndFlux_c,bndFlux_v)}
        else:
            boundary = {(1,2): Psi, (3,4): (bndFlux_c,bndFlux_v)}
    return Model

# 1) Dirichlet boundary conditions with exact solution
Psi = test(dirichletTest,None,None)
# 2) now use discrete solution from first step as boundary condition
testPsi = test(dirichletTest,Psi, show=False, tol=0)
# 3) now use Neumann conditions based on exact solution
testPsi = test(lambda p: neumannTest(p,exact=True),Psi, show=False, tol=0)
testPsi = test(lambda p: neumannTest(p,exact=False),Psi, show=False, tol=0.001)
