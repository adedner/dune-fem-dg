from ._operators import *

from .model2ufl import model2ufl, model2dgufl
from .boundary import BndValue, BndFlux_v, BndFlux_c

registry = {}

registry["scheme"] = {
        "rungekutta": rungeKuttaSolver
        }

def _cite_dune_module_as():
    return """
@article{dunefemdg:21,
  title={{Extendible and Efficient Python Framework for Solving Evolution Equations
          with Stabilized Discontinuous Galerkin Methods}},
  author={Dedner, A. and Kl{\\"o}fkorn, R.},
  year={2021},
  journal = {{ Commun. Appl. Math. Comput.}},
  doi={10.1007/s42967-021-00134-5}
}
"""
