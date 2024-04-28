import os

# set number of threads to be used for thread parallel version
os.environ['OMP_NUM_THREADS'] = '4'

from dune.fem import parameter
from dune.femdg.testing import run

# from scalar import shockTransport as problem
# from scalar import sinProblem as problem
# from scalar import sinTransportProblem as problem
# from scalar import sinAdvDiffProblem as problem
from scalar import pulse as problem
# from scalar import diffusivePulse as problem

parameter.append({"fem.verboserank": 0})

parameters = {"fem.ode.odesolver": "EX",   # EX, IM, IMEX
              "fem.ode.order": 3,
              "fem.ode.verbose": "cfl",      # none, cfl, full
              "fem.timeprovider.factor": 0.45,
              "dgadvectionflux.method": "LLF",
              "dgdiffusionflux.method": "CDG2",      # CDG2, CDG, BR2, IP, NIPG, BO
              "dgdiffusionflux.theoryparameters": 1, # scaling with theory parameters
              "dgdiffusionflux.penalty": 0,
              "dgdiffusionflux.liftfactor": 1}

parameters['fem.ode.odesolver'] = 'EX'
Model = problem()
Model.endTime = 0.1
uh,errorEx = run(Model,
        startLevel=0, polOrder=2, limiter='scaling',
        primitive=None, saveStep=0.01, threading=True, grid="alucube", subsamp=0, space="dgonb",
        dt=None,
        parameters=parameters)
