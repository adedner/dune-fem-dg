import math
from ufl import *
from dune.grid import structuredGrid
from dune.fem.space import dgonb
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper

# Basic model for hyperbolic conservation law
class Model:
    gamma = 1.4
    # helper function
    def toPrim(U):
        v = as_vector( [U[i]/U[0] for i in range(1,3)] )
        kin = dot(v,v) * U[0] / 2
        pressure = (Model.gamma-1)*(U[3]-kin)
        return U[0], v, pressure

    # interface methods for model
    def F_c(t,x,U):
        rho, v, p = Model.toPrim(U)
        return as_matrix( [
                  [rho*v[0], rho*v[1]],
                  [rho*v[0]*v[0] + p, rho*v[0]*v[1]],
                  [rho*v[0]*v[1], rho*v[1]*v[1] + p],
                  [(U[3]+p)*v[0], (U[3]+p)*v[1]] ] )
    # simple 'outflow' boundary conditions on all boundaries
    boundary = {range(1,5): lambda t,x,U: U}

    # interface method needed for LLF and time step control
    def maxLambda(t,x,U,n):
        rho, v, p = Model.toPrim(U)
        return abs(dot(v,n)) + sqrt(Model.gamma*p/rho)

########################################################################

# Part 1: basic set up and time loop - no limiter and fixed time step
#         using constant initial data
gridView = structuredGrid([-1,-1],[1,1],[40,40])
space = dgonb( gridView, order=3, dimRange=4)
operator = femDGOperator(Model, space, limiter=None)
print(operator.__module__)
stepper  = femdgStepper(order=3, operator=operator)
u_h = space.interpolate([1.4,0,0,1], name='u_h')
t  = 0
# don't show vtk in paper
vtk = gridView.sequencedVTK("paperA", pointdata=[u_h])
while t < 0.1:
    operator.setTime(t)
    t += stepper(u_h, dt=0.001)
vtk()

# Part 2: add limiter and methods for troubled cell indicator
def velocity(t,x,U):
    _, v ,_ = Model.toPrim(U)
    return v
def physical(U):
    rho, _, p = Model.toPrim(U)
    return conditional( rho>1e-8, conditional( p>1e-8 , 1, 0 ), 0 )
def jump(U,V):
    _,_, pL = Model.toPrim(U)
    _,_, pR = Model.toPrim(V)
    return (pL - pR)/(0.5*(pL + pR))
# don't show the following lines just explain that the above method is added to the Model
Model.velocity = velocity
Model.physical = physical
Model.jump     = jump

operator = femDGOperator(Model, space, limiter="MinMod")
print(operator.__module__)
stepper  = femdgStepper(order=3, operator=operator, cfl=0.2)
x = SpatialCoordinate(space)
u_h.interpolate( conditional(dot(x,x)<0.1,as_vector([1,0,0,2.5]),
                                          as_vector([0.125,0,0,0.25])) )
operator.applyLimiter(u_h)
t  = 0
vtk = gridView.sequencedVTK("paperB", pointdata=[u_h])
while t < 0.4:
    vtk()
    assert not math.isnan( u_h.scalarProductDofs( u_h ) )
    operator.setTime(t)
    t += stepper(u_h)
    print(t)
vtk()
