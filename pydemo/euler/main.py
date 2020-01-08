import mpi4py.rc
mpi4py.rc.threaded = True
from dune.fem import parameter
from dune.femdg.testing import run

from euler import sod as problem
# from euler import vortex as problem
# from euler import leVeque as problem
# from euler import radialSod3 as problem

dim = 2
gamma = 1.4

parameter.append({"fem.verboserank": -1})

primitive=lambda Model,uh: {"pressure": Model.toPrim(uh)[2]}
parameters = {"fem.ode.odesolver": "EX",
              "fem.timeprovider.factor": 0.25,
              "fem.ode.order": 3,
              "dgadvectionflux.method": "EULER-HLLC",
              "femdg.limiter.admissiblefunctions": 1,
              "femdg.limiter.tolerance": 1,
              "femdg.limiter.epsilon": 1e-8}
#-----------------
# femdg.limiter.tolerance: 1 (tolerance for shock indicator)
# femdg.limiter.epsilon: 1e-8 (epsilon to avoid rounding errors)
# femdg.limiter.admissiblefunctions:
#    0 = only dg solution | 1 = only reconstruction | 2 = both
#-----------------


run(problem(),
    startLevel=0, polOrder=2, limiter="default",
    primitive=primitive, saveStep=0.16, subsamp=2,
    dt=None,threading=True,grid="alucube", space="dgonb",
    parameters=parameters)
