from dune.fem import parameter
from dune.femdg.testing import run

from euler import sod as problem
# from euler import vortex as problem
# from euler import leVeque as problem
# from euler import radialSod3 as problem

dim = 2
gamma = 1.4

parameter.append({"fem.verboserank": -1})
parameter.append("parameter")

primitive=lambda Model,uh: {"freesurface":Model.toPrim(uh)[0]}

run(*problem(dim,gamma),
        startLevel=0, polOrder=2, limiter="default",
        primitive=primitive, saveStep=0.01, subsamp=2)
