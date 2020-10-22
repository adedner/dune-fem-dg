import os
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('problem', type=int)
parser.add_argument
parser.parse_args()
try:
    args = parser.parse_args()
except SystemExit:
    sys.exit(0)
problem = args.problem

from dune.alugrid import aluCubeGrid as aluGrid
from dune.grid import reader, cartesianDomain
from dune.fem.view import adaptiveLeafGridView as view
from dune.polygongrid import polygonGrid

import dune.fem
# dune.fem.parameter.append({"fem.verboserank":0})
# dune.fem.parameter.append({"fem.parallel.numberofthreads":8})

from advection_model import model
from evolve import evolve

Model = model(2)
order = 4
domain = (reader.dgf, "unitsquare.dgf")
os.makedirs("advection", exist_ok=True)

##################

if problem == 1:
    gridView = view( aluGrid( domain, dimgrid=2 ) )
    gridView.hierarchicalGrid.globalRefine(3)
    evolve(gridView, order, Model, "advection/none",
           space="onb", limiter=None,
           maxLevel=-1, threading=False, codegen=True, outputs=5)

##################

if problem == 2:
    gridView = view( aluGrid( domain, dimgrid=2 ) )
    gridView.hierarchicalGrid.globalRefine(3)
    evolve(gridView, order, Model, "advection/minmod",
           space="onb", maxLevel=-1, threading=False, codegen=True, outputs=5)

if problem == 3:
    from ufl import conditional
    delattr(Model,"jump")
    delattr(Model,"velocity")
    Model.physical = lambda t,x,U: (
              conditional( U[0]>-1e-8, 1.0, 0.0 ) *
              conditional( U[0]<1.0+1e-8, 1.0, 0.0 ) )
    gridView = view( aluGrid( domain, dimgrid=2 ) )
    gridView.hierarchicalGrid.globalRefine(3)
    evolve(gridView, order, Model, "advection/physical",
           space="lobatto", limiter="lp",
           maxLevel=-1, threading=False, codegen=True, outputs=5)

if problem == 4:
    from ufl import conditional
    delattr(Model,"jump")
    delattr(Model,"velocity")
    Model.lowerBound = [0]
    Model.upperBound = [1]
    gridView = view( aluGrid( domain, dimgrid=2 ) )
    gridView.hierarchicalGrid.globalRefine(3)
    evolve(gridView, order, Model, "advection/scaling",
           #space="lobatto", limiter="scaling",
           space="legendre", limiter="scaling",
           maxLevel=-1, threading=False, codegen=True, outputs=5)

##################

if problem == 5:
    gridView = view( aluGrid( domain, dimgrid=2 ) )
    evolve(gridView, order, Model, "advection/adapt_none",
           space="lobatto",
           limiter=None, maxLevel=5, codegen=True)

##################

if problem == 6:
    gridView = view( aluGrid( domain, dimgrid=2 ) )
    evolve(gridView, order, Model, "advection/adapt_limit",
           maxLevel=5, codegen=True)

##################

# domain = (reader.dgf, "triangle.dgf")
# gridView = view( polygonGrid( domain, dualGrid=True ) )
# parameters["femdg.limiter.admissiblefunctions"]=3
# evolve(gridView, order, Model, "advection/poly",  limiter="MinMod",
#              maxLevel=-1)
