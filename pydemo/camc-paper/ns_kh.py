import os, sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('level', type=int)
parser.add_argument('--space', type=str, default="onb")
parser.add_argument('--stepper', type=str, default="femdg")
parser.add_argument('--order', type=int, default=4)
parser.add_argument('--grid',  type=str, default="cube")
parser.add_argument('--out',   type=int, default=100)
parser.add_argument('--mu',   type=float, default=0.001)
parser.add_argument
parser.parse_args()
try:
    args = parser.parse_args()
except SystemExit:
    sys.exit(0)
space = args.space
stepper = args.stepper
level = args.level
dim   = 2
grid  = args.grid
order = args.order
out   = args.out
mu    = args.mu

path = "euler_tmp/"
os.makedirs(path, exist_ok=True)
path = path+grid+str(dim)+str(abs(level))+"_"+space+"_"+str(mu).replace(".","")+"_"
if stepper != "femdg": path = path+stepper

from dune.grid import reader
from ns_model import model
Model = model(problem="KH",dim=2, mu=mu)

from dune.fem import parameter
parameter.append({"fem.verboserank":0})

if grid == "cube":
    from dune.alugrid import aluCubeGrid as grid
elif grid == "simplex":
    from dune.alugrid import aluSimplexGrid as grid
elif grid == "naffine":
    from dune.alugrid import aluCubeGrid as grid
    Model.domain = (reader.dgf, "shockvortex_naffine.dgf")
elif grid == "poly":
    from dune.polygongrid import polygonGrid as gridPoly
    Model.domain = (reader.dgf, "shockvortex_poly.dgf")
    Model.boundary = Model.polyBoundary
    grid = lambda domain, dimgrid: gridPoly( Model.domain, dualGrid=True )

from dune.fem.view import adaptiveLeafGridView as view
from evolve import evolve

gridView = view( grid( Model.domain, dimgrid=dim ) )

parameters = {"fem.ode.verbose":"full",
              "fem.solver.newton.verbose":True,
              "fem.solver.newton.linear.verbose":True,
              "fem.solver.newton.linear.gmres.restart":100}


if level <= 0:
    if level<0:
        gridView.hierarchicalGrid.globalRefine(-level)
    evolve(gridView, order, Model, path,
           limiter=None,
           maxLevel=-1, outputs=out,
           parameters=parameters,
           space=space, stepper=stepper)
else:
    # gridView.hierarchicalGrid.globalRefine(level//2)
    evolve(gridView, order, Model, path+"_adapt",
           limiter=None,
           maxLevel=level, outputs=out,
           parameters=parameters,
           space=space, stepper=stepper)
