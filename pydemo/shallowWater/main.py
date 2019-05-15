from dune.fem import parameter
from dune.femdg.testing import run

from shallowwater import leVeque as problem

parameter.append({"fem.verboserank": -1})
parameter.append("parameter")

primitive=lambda Model,uh: {"freesurface":Model.toPrim(uh)[0]}
parameters = {"fem.ode.odesolver": "EX",
              "fem.timeprovider.factor": 0.45/2.,
              "dgadvectionflux.method": "LLF",
              "femdg.limiter.limiteps": 1,
              "femdg.limiter.admissiblefunctions": 1,
              "femdg.limiter.tolerance": 1}
run(*problem(1),
        startLevel=0, polOrder=2, limiter="default",
        primitive=primitive, saveStep=0.1, subsamp=0,
        dt=None,threading=True, grid="alucube",
        parameters=parameters)
