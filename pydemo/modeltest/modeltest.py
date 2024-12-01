import dune.femdg

import numpy as np
from matplotlib import pyplot
from ufl import *
from dune.ufl import DirichletBC, Constant
from dune.grid import structuredGrid
from dune.fem.function import gridFunction
from dune.fem.plotting import plotPointData as plot
from dune.fem.space import lagrange, dgonb
from dune.fem.scheme import galerkin
from dune.femdg import femDGUFL
from dune.femdg.rk import femdgStepper

from dune.fem import threading
threading.use = max(4,threading.max) # use at most 4 threads

# Given grid function u_n and
# L[U] = L_e[U] + L_i[U] = (-div(F_c) - S_e)) + (-div(F_v) - S_i)
# implement model for operator
# T[U] = (U - u_n) + tau L[U]
# So with F_v = grad(u), F_c=S_e=0:
# T[U] = U - u_n - laplace(U) - S_i
# so we modify S_i <= S_i - U and S_e <= S_e + u_n

def timeModel(Model, old):
    class TimeModel(Model):
        tau = dune.ufl.Constant(0,"tau")
        if hasattr(Model,"S_i"):
            def S_i(t,x,U,DU): # or S_e for a non stiff source
                return TimeModel.tau * Model.S_i(t,x,U,DU) - U
        else:
            def S_i(t,x,U,DU): # or S_e for a non stiff source
                return -U
        if hasattr(Model,"S_e"):
            def S_e(t,x,U,DU): # or S_e for a non stiff source
                return TimeModel.tau * Model.S_e(t,x,old,grad(old)) + old
        else:
            def S_e(t,x,U,DU): # or S_e for a non stiff source
                return old
        if hasattr(Model,"F_c"):
            def F_c(t,x,U):
                return TimeModel.tau * Model.F_c(t,x,old)
        if hasattr(Model,"F_v"):
            def F_v(t,x,U,DU):
                return TimeModel.tau * Model.F_v(t,x,U,DU)
    return TimeModel

################################################################

class VelocityModel:
    def S_e(t,x,U,DU):
        return conditional(x[1]<0,as_vector([3]),as_vector([-3]))
    def F_v(t,x,U,DU):
        return DU
    boundary = {range(1,5): as_vector([0.])}

gridView = structuredGrid([0,-np.pi],[np.pi,np.pi],[30,60])
space = lagrange(gridView,dimRange=1)
psi = space.interpolate(0,name="streamFunction")
scheme = galerkin( femDGUFL(VelocityModel, space) )
scheme.solve(target=psi)
velocity = as_vector([-psi[0].dx(1),psi[0].dx(0)])
# gridFunction(velocity).plot(gridLines=None, vectors=[0,1], block=False)
# pyplot.show()

#############################################################

class ChemicalModel:
    def S_e(t,x,U,DU):
        P1 = as_vector([0.2*pi,-0.8*pi]) # midpoint of first source
        P2 = as_vector([0.2*pi, 0.8*pi]) # midpoint of second source
        f1 = conditional(dot(x-P1,x-P1) < 0.2, 1, 0)
        f2 = conditional(dot(x-P2,x-P2) < 0.2, 1, 0)
        f  = conditional(t<5, as_vector([f1,f2,0]), as_vector([0,0,0]))
        r = 10*as_vector([U[0]*U[1], U[0]*U[1], -2*U[0]*U[1]])
        return (f - r) # issue here - other code has this as (f-r)
    def F_c(t,x,U):
        return as_matrix([ [*(velocity*u)] for u in U ])
    def F_v(t,x,U,DU):
        return 0.02*DU
    boundary = {range(1,5): as_vector([0,0,0])}

space = lagrange(gridView,dimRange=3)
c_old = space.interpolate(as_vector([0,0,0]),name="concentrations")
c_new = space.interpolate(as_vector([0,0,0]),name="concentrations")
Model = timeModel(ChemicalModel, c_old)
scheme = galerkin( femDGUFL(Model, space) )

###################################

tau = 0.01
saveInterval = 100*tau
nextSaveTime = saveInterval
Model.tau.value = tau
endTime = 10
t = 0
while t<endTime:
    c_old.assign(c_new)
    scheme.model.time = t
    scheme.solve(target=c_new)
    t += tau
    if t > nextSaveTime:
        print("time=",t,flush=True)
        # c_new.plot()
        nextSaveTime += saveInterval

from dune.fem.plotting import plotComponents
from matplotlib import ticker
# plotComponents(c_new, gridLines=None, level=1,
#                colorbar={"orientation":"horizontal", "ticks":ticker.MaxNLocator(nbins=4)})
