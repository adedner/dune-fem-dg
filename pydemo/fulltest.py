from dune.fem import parameter
from dune.femdg.testing import run

import scalar, shallowWater, euler

parameter.append({"fem.verboserank": -1})
parameter.append("parameter")

for p in shallowWater.problems + euler.problems:
    run(*p(), startLevel=0, polOrder=2, limiter="default",
              primitive=None, saveStep=None, subsamp=0)
