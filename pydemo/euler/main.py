import mpi4py.rc
mpi4py.rc.threaded = True
from dune.fem import parameter
from dune.femdg.testing import run

# from euler import sod as problem
from euler import vortex as problem
# from euler import leVeque as problem
# from euler import radialSod3 as problem

dim = 2
gamma = 1.4

parameter.append({"fem.verboserank": -1})
parameter.append({"fem.parallel.numberofthreads": 4})
parameter.append("parameter")

primitive=lambda Model,uh: {"pressure":Model.toPrim(uh)[2]}

run(*problem(),
        startLevel=0, polOrder=2, limiter="default",
        primitive=primitive, saveStep=0.1, subsamp=2,
        dt=None,threading=True,grid="yasp")
