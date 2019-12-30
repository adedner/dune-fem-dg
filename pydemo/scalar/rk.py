from dune.fem import parameter
from dune.femdg.rk import run, femdgStepper,\
        Heun,explSSP3,ExplSSP4_10

# from scalar import shockTransport as problem
# from scalar import sinProblem as problem
# from scalar import sinTransportProblem as problem
# from scalar import sinAdvDiffProblem as problem
from scalar import pulse as problem
# from scalar import diffusivePulse as problem

parameter.append({"fem.verboserank": 0})
parameters = {"dgadvectionflux.method": "LLF",
              "fem.solver.gmres.restart": 50,
              "dgdiffusionflux.method": "CDG2",      # CDG2, CDG, BR2, IP, NIPG, BO
              "dgdiffusionflux.theoryparameters": 1, # scaling with theory parameters
              "dgdiffusionflux.penalty": 0,
              "dgdiffusionflux.liftfactor": 1,
              "fem.ode.odesolver": "EX",   # EX, IM, IMEX
              "fem.ode.order": 3,
              "fem.ode.verbose": "cfl",      # none, cfl, full
              "fem.ode.cflMax": 100,
              "fem.ode.cflincrease": 1.25,
              "fem.ode.miniterations": 35,
              "fem.ode.maxiterations": 100}

startLevel=1
uh,error = run(problem(), femdgStepper(parameters),
        startLevel=startLevel, polOrder=2, limiter=None,
        primitive=None, saveStep=0.01, subsamp=0,
        dt=None, parameters=parameters)
uh,errorHeun = run(problem(), Heun,
        startLevel=startLevel, polOrder=2, limiter=None,
        primitive=None, saveStep=0.01, subsamp=0,
        dt=None, parameters=parameters)
uh,errorSSP3 = run(problem(), explSSP3(4),
        startLevel=startLevel, polOrder=2, limiter=None,
        primitive=None, saveStep=0.01, subsamp=0,
        dt=None, parameters=parameters)
uh,errorSSP4 = run(problem(), ExplSSP4_10,
        startLevel=startLevel, polOrder=2, limiter=None,
        primitive=None, saveStep=0.01, subsamp=0,
        dt=None, parameters=parameters)

print(error,errorHeun,errorSSP3,errorSSP4)
