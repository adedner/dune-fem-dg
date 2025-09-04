import inspect

from ufl import replace, zero
from ufl.algorithms.ad import expand_derivatives
from ufl.core.expr import Expr

class BoundaryCondition:
    """ base class for boundary conditions in Model.boundary
    """
    def __init__(self, value):
        self.value = value
    def __call__(self, *args, **kwds):
        return self.value(*args, **kwds)


class BndValue(BoundaryCondition):
    """ class for Dirichlet type boundary condition

    value can be an expression, or a callable with
    - one argument (x)
    - two arguments (t,x)
    - three arguments (t,x,u)
    """
    def __init__(self, value):
        if isinstance(value, Expr):
            super().__init__(lambda t, x, u: value)
        else:
            num_args = len(inspect.signature(value).parameters)
            if num_args == 1:
                super().__init__(lambda t, x, u: value(x))
            elif num_args == 2:
                super().__init__(lambda t, x, u: value(t, x))
            elif num_args == 3:
                super().__init__(value)
            else:
                raise ValueError(f"Boundary has {num_args} arguments.")


class BndFlux_v(BoundaryCondition):
    """ class for viscous flux boundary, 'value' should be callable
    with signmature (t,x,u,Du,n)
    """
    pass


class BndFlux_c(BoundaryCondition):
    """ class for convective flux boundary, 'value' should be callable
    with signmature (t,x,u,n)
    """
    pass


def classify_boundary(boundary_dict):
    """ utility method that splits a boundary dictionary into three parts
    for Dirichlet,convective and diffusive fluxes
    """
    boundary_flux_cs = {}  # Fluxes for the advection term
    boundary_flux_vs = {}  # Fluxes for the diffusion term
    boundary_values = {}   # Boundary values for Dirichlet

    for k, f in boundary_dict.items():

        if isinstance(k, (Expr, str)):
            ids = [k]
        elif callable(k):
            ids = [k]
        else:
            try:
                ids = []
                for kk in k:
                    ids += [kk]
            except TypeError:
                ids = [k]

        if isinstance(f, (tuple, list)):
            assert len(f) == 2, "too many boundary fluxes provided"
            if isinstance(f[0], BndFlux_v) and isinstance(f[1], BndFlux_c):
                boundary_flux_vs.update([(kk, f[0]) for kk in ids])
                boundary_flux_cs.update([(kk, f[1]) for kk in ids])

            elif isinstance(f[0], BndFlux_c) and isinstance(f[1], BndFlux_v):
                boundary_flux_vs.update([(kk, f[1]) for kk in ids])
                boundary_flux_cs.update([(kk, f[0]) for kk in ids])

            elif isinstance(f[1], BndFlux_c) and isinstance(f[0], BndFlux_v):
                boundary_flux_vs.update([(kk, f[0]) for kk in ids])
                boundary_flux_cs.update([(kk, f[1]) for kk in ids])

            else:
                raise ValueError("Need AFlux and DFlux")

        elif isinstance(f, BndFlux_v):
            boundary_flux_vs.update([(kk, f) for kk in ids])

        elif isinstance(f, BndFlux_c):
            boundary_flux_cs.update([(kk, f) for kk in ids])

        elif isinstance(f, BndValue):
            boundary_values.update([(kk, f) for kk in ids])

        else:
            raise NotImplementedError(f"unknown boundary type {k} : {f}")

    return boundary_flux_cs, boundary_flux_vs, boundary_values


def splitBoundary(Model, override_boundary_dict=None):
    """ take a Model and split boundary dictionary testing and requirements are met
    """
    if override_boundary_dict is not None:
        boundary_dict = override_boundary_dict
    else:
        boundary_dict = Model.boundary

    hasFlux_c = hasattr(Model, "F_c")
    hasFlux_v = hasattr(Model, "F_v")

    boundary_flux_cs, boundary_flux_vs, boundary_values = classify_boundary(
        boundary_dict
    )

    if hasFlux_c and hasFlux_v:
        assert len(boundary_flux_cs) == len(
            boundary_flux_vs
        ), "two bulk fluxes given, but one boundary fluxes provided"

    if not hasFlux_c:
        assert len(boundary_flux_cs) == 0, "No bulk Advection, but boundary flux given"

    if not hasFlux_v:
        assert len(boundary_flux_vs) == 0, "No bulk diffusion, but boundary flux given"

    assert boundary_values.keys().isdisjoint(boundary_flux_cs)
    assert boundary_values.keys().isdisjoint(boundary_flux_vs)

    return boundary_flux_cs, boundary_flux_vs, boundary_values
