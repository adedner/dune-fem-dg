from dune.fem import parameter
from dune.femdg.testing import run

from shallowwater import leVeque as problem

parameter.append({"fem.verboserank": -1})
parameter.append("parameter")

primitive=lambda Model,uh: {"freesurface":Model.toPrim(uh)[0]}

run(*problem(1),
        startLevel=0, polOrder=2, limiter="default",
        primitive=primitve, saveStep=0.1, subsamp=0)
