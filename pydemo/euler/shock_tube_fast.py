import time
from dune.grid import structuredGrid, cartesianDomain
from dune.fem import parameter
import dune.create as create
from dune.models.elliptic.formfiles import loadModels
from llf import NumFlux
from dune.femdg import femDGOperator, rungeKuttaSolver
from ufl import *

gamma = 1.4
dim = 2

from euler import sod as problem
#from euler import radialSod3 as problem

Model, initial, domain, endTime, name, exact = problem(dim,gamma)

parameter.append({"fem.verboserank": -1})
parameter.append({"fem.timeprovider.factor": 0.25})
parameter.append("parameter")

parameters = {"fem.ode.odesolver": "EX",
              "fem.timeprovider.factor": "0.35",
              "dgadvectionflux.method": "EULER-HLLC",
              "femdg.limiter.admissiblefunctions": 1,
              "femdg.limiter.tolerance": 1e-8}

x0,x1,N = domain
grid = structuredGrid(x0,x1,N)
# grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
dimR     = grid.dimension + 2
t        = 0
dt       = 1e-3
count    = 0
saveStep = 0.15
saveTime = saveStep

def initialize(space):
    return space.interpolate(initial, name='u_h')
    if space.order == 0:
        return space.interpolate(initial, name='u_h')
    else:
        lagOrder = 1 # space.order
        spacelag = create.space("lagrange", space.grid, order=lagOrder, dimRange=space.dimRange)
        u_h = spacelag.interpolate(initial, name='tmp')
        return space.interpolate(u_h, name='u_h')

def useODESolver(polOrder=2, limiter='default'):
    global count, t, dt, saveTime
    polOrder = polOrder
    if False:
        # needs change in dune/fem-dg/operator/dg/passtraits.hh
        space = create.space("dglegendre", grid, order=polOrder, dimRange=dimR, hierarchical=False)
    else:
        space = create.space("dgonb", grid, order=polOrder, dimRange=dimR)
    u_h = initialize(space)
    # rho, v, p = Model.toPrim(u_h)
    operator = femDGOperator(Model, space, limiter=limiter, threading=True, parameters=parameters )
    ode = rungeKuttaSolver( operator, parameters=parameters )

    operator.applyLimiter( u_h )
    print("number of elements: ",grid.size(0),flush=True)
    grid.writeVTK(name,
        pointdata=[u_h],
        # celldata={"density":rho, "pressure":p}, # bug: density not shown correctly
        #celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
        #cellvector={"velocity":v},
        number=count, subsampling=2)
    start = time.time()
    tcount = 0
    while t < endTime:
        ode.solve(u_h)
        # operator.applyLimiter( u_h );
        dt = ode.deltaT()
        print('dt = ',dt)
        t += dt
        tcount += 1
        if tcount%100 == 0:
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
        if t > saveTime:
            count += 1
            grid.writeVTK(name,
                pointdata=[u_h],
                #celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
                #cellvector={"velocity":v},
                number=count, subsampling=2)
            saveTime += saveStep
    print("time loop:",time.time()-start)
    print("number of time steps ", tcount)
    grid.writeVTK(name,
        pointdata=[u_h],
        #celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
        #cellvector={"velocity":v},
        number=count, subsampling=2)

scheme = 2

if scheme == 0:
    # grid = structuredGrid(x0,x1,N)
    grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    #grid = create.grid("ALUCube", cartesianDomain(x0,x1,N))
    grid.hierarchicalGrid.globalRefine(1)
    # grid = create.view("adaptive", grid)
    useODESolver(2,'default')      # third order with limiter
elif scheme == 1:
    #N = [n*4 for n in N]
    #grid = structuredGrid(x0,x1,N)
    grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    #grid = create.grid("ALUCube", cartesianDomain(x0,x1,N))
    grid.hierarchicalGrid.globalRefine(1)
    # grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    useODESolver(0,None)           # FV scheme
elif scheme == 2:
    grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    #grid = create.grid("ALUCube", cartesianDomain(x0,x1,N))
    grid.hierarchicalGrid.globalRefine(1)
    #N = [n*6 for n in N]
    #grid = structuredGrid(x0,x1,N)
    # grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    useODESolver(0,'default')      # FV scheme with limiter
