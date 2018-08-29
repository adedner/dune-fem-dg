from dune.fem import parameter
from dune.femdg.testing import run

# from scalar import shockTransport as problem
# from scalar import sinProblem as problem
from scalar import sinTransportProblem as problem
# from scalar import pulse as problem
# from scalar import diffusivePulse as problem

parameter.append({"fem.verboserank": 0})
parameter.append("parameter")

run(*problem(),
        startLevel=0, polOrder=2, limiter=None,
        primitive=None, saveStep=0.1, subsamp=0)
