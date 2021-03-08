from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-n','--num-threads', type=int,
        help="""
        Number of threads to use when running the test""")
parser.add_argument('-s','--scheme', type=int, default=0,
        help="""
        Scheme: 0 == DG | 1 == FV0 | 2 == FV1""")
parser.add_argument('-o','--order', type=int, default=2,
        help="""
        Polynomial order of DG scheme (default = 2)""")
parser.add_argument
parser.parse_args()
try:
    args = parser.parse_args()
except SystemExit:
    sys.exit(0)

nThreads = args.num_threads
scheme = args.scheme
dgOrder = args.order
print(f"Running with scheme {scheme} with {nThreads} threads!")

import os
# set number of threads to be used for thread parallel version
os.environ['OMP_NUM_THREADS'] = str(nThreads)

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

Model = problem(dim,gamma)

parameter.append({"fem.verboserank": 0})
parameter.append({"fem.timeprovider.factor": 0.35})
parameter.append({"fem.threads.verbose" : True })
parameter.append({"fem.adaptation.method" : "none"})

#parameter.append("parameter")

parameters = {"fem.ode.odesolver": "EX",
              "fem.ode.order" : 3,
              "fem.timeprovider.factor": 0.35,
#              "femdg.limiter.admissiblefunctions" : 1,
              "femdg.nonblockingcomm" : True
 #            "dgadvectionflux.method": "EULER-LLF",
        }

parameter.append(parameters)

x0,x1,N = Model.domain
# grid = structuredGrid(x0,x1,N)
# grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
dimR = Model.dimRange
t        = 0
count    = 0
saveStep = 0.15
saveTime = saveStep

def initialize(space):
    return space.interpolate(Model.initial, name='u_h')
    #if space.order == 0:
    #    return space.interpolate(initial, name='u_h')
    #else:
    #    lagOrder = 1 # space.order
    #    spacelag = create.space("lagrange", space.grid, order=lagOrder, dimRange=space.dimRange)
    #    u_h = spacelag.interpolate(initial, name='tmp')
    #    return space.interpolate(u_h, name='u_h')

def useODESolver(polOrder=2, limiter='default', codegen=True, spc='onb'):
    global count, t, dt, saveTime
    polOrder = polOrder
    if spc == 'lobatto':
        from dune.fem.space import dglagrange
        space = dglagrange( grid, order=polOrder, dimRange=dimR, pointType='lobatto', codegen=False)
        #space = create.space("dglegendre", grid, order=polOrder, dimRange=dimR, hierarchical=False)
    else:
        from dune.fem.space import dgonb
        space = dgonb(grid, order=polOrder, dimRange=dimR, codegen=False)
    u_h = initialize(space)
    # rho, v, p = Model.toPrim(u_h)
    operator = femDGOperator(Model, space, limiter=limiter, threading=False )
    ode = rungeKuttaSolver( operator )

    operator.applyLimiter( u_h )
    print("number of elements: ",grid.size(0),flush=True)
    grid.writeVTK(Model.name,
        pointdata=[u_h],
        # celldata={"density":rho, "pressure":p}, # bug: density not shown correctly
        #celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
        #cellvector={"velocity":v},
        number=count, subsampling=2)
    start = time.time()
    tcount = 0

    #dt = 0.0001325
    while t < Model.endTime:
        #ode.setTimeStepSize( dt )
        ode.solve(u_h)
        dt = ode.deltaT()
        t += dt
        tcount += 1
        if tcount%100 == 0:
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
        if t > saveTime:
            count += 1
            grid.writeVTK(Model.name,
                pointdata=[u_h],
                #celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
                #cellvector={"velocity":v},
                number=count, subsampling=2)
            saveTime += saveStep
        if t + dt > Model.endTime:
            print('[',tcount,']','dt = ', dt, 'time = ',t, 'count = ',count, flush=True )
            break

    print("time loop:",time.time()-start)
    print("number of time steps ", tcount)
    grid.writeVTK(Model.name,
        pointdata=[u_h],
        #celldata={"pressure":p, "maxLambda":Model.maxLambda(0,0,u_h,as_vector([1,0]))},
        #cellvector={"velocity":v},
        number=count, subsampling=2)

scheme = 0

if scheme == 0:
    # grid = structuredGrid(x0,x1,N)
    #grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    grid = create.grid("ALUCube", cartesianDomain(x0,x1,N))
    grid.hierarchicalGrid.globalRefine(3)
    # grid = create.view("adaptive", grid)
    useODESolver(dgOrder,'default',spc='onb')      # third order with limiter
elif scheme == 1:
    #N = [n*4 for n in N]
    #grid = structuredGrid(x0,x1,N)
    grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    #grid = create.grid("ALUCube", cartesianDomain(x0,x1,N))
    grid.hierarchicalGrid.globalRefine(1)
    # grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    useODESolver(0,None)           # FV scheme
elif scheme == 2:
    #grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    grid = create.grid("ALUCube", cartesianDomain(x0,x1,N))
    grid.hierarchicalGrid.globalRefine(1)
    #N = [n*6 for n in N]
    #grid = structuredGrid(x0,x1,N)
    # grid = create.grid("ALUSimplex", cartesianDomain(x0,x1,N))
    useODESolver(0,'default')      # FV scheme with limiter

parameter.write("param.log")
