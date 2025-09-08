from ._operators import *

from .model2ufl import model2ufl
from .boundary import BndValue, BndFlux_v, BndFlux_c

registry = {}

registry["scheme"] = {
        "rungekutta": rungeKuttaSolver
        }
