import os

# set number of threads to be used for thread parallel version
os.environ['OMP_NUM_THREADS'] = '4'

from dune.fem import parameter
from dune.femdg.testing import run

# from euler import constant as problem
from euler import sod as problem
# from euler import vortex as problem
# from euler import leVeque as problem
# from euler import radialSod3 as problem

dim = 2
gamma = 1.4

parameter.append({"fem.verboserank": 0})

primitive=lambda Model,uh: {"pressure": Model.toPrim(uh)[2]}
parameters = {"fem.ode.odesolver": "EX",
              "fem.ode.order": 3,
              "femdg.limiter.tolerance":1 }

#-----------------
# "dgadvectionflux.method": "EULER-HLLC", "EULER-HLL", "LLF"
# default value is 'LLF'
#-----------------
# femdg.limiter.tolerance: 1 (tolerance for shock indicator)
# femdg.limiter.epsilon: 1e-8 (epsilon to avoid rounding errors)
# femdg.limiter.admissiblefunctions:
#    0 = only dg solution | 1 = only reconstruction | 2 = both
#-----------------

Model = problem()
#Model.endTime = 0.1501
# Model.exact = None

run(Model,
    startLevel=0, polOrder=2, limiter="default",
    primitive=None, saveStep=0.16, subsamp=0,
    dt=None,threading=True,grid="yasp",
    space="dgonb",
    #space="dglagrange",
    #space=("dglagrange","lobatto"),
    # space=("dglagrange","gauss"),
    parameters=parameters)

# print(str(parameter))
