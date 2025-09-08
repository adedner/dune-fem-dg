import inspect

import ufl
from ufl.core.expr import Expr
from dune.fem.deprecated import deprecated

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
    - four arguments (t,x,u,n)
    """
    def __init__(self, value):
        if isinstance(value, Expr):
            super().__init__(lambda t,x,u,n: value)
        else:
            num_args = len(inspect.signature(value).parameters)
            if num_args == 1:
                super().__init__(lambda t,x,u,n: value(x))
            elif num_args == 2:
                super().__init__(lambda t,x,u,n: value(t, x))
            elif num_args == 3:
                super().__init__(lambda t,x,u,n: value(t, x, u))
            elif num_args == 4:
                super().__init__(lambda t,x,u,n: value(t, x, u, n))
            elif num_args == 5: # this version is used in old setup
                # and should already produce a deprecation warning
                super().__init__(lambda t,x,u,n: value(t, x, u, n, n))
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


def classify_boundary(boundary_dict, hasAdvFlux, hasDiffFlux):
    """ utility method that splits a boundary dictionary into three parts
    for Dirichlet,convective and diffusive fluxes
    """
    boundary_flux_cs = {}  # Fluxes for the advection term
    boundary_flux_vs = {}  # Fluxes for the diffusion term
    boundary_values = {}   # Boundary values for Dirichlet

    for k, f in boundary_dict.items():

        # collect all ids
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

        needOld = False
        if isinstance(f, (tuple, list)):
            assert len(f) == 2, "too many boundary fluxes provided"
            if isinstance(f[0], BndFlux_v) and isinstance(f[1], BndFlux_c):
                boundary_flux_vs.update([(kk, f[0]) for kk in ids])
                boundary_flux_cs.update([(kk, f[1]) for kk in ids])

            elif isinstance(f[0], BndFlux_c) and isinstance(f[1], BndFlux_v):
                boundary_flux_vs.update([(kk, f[1]) for kk in ids])
                boundary_flux_cs.update([(kk, f[0]) for kk in ids])

            else:
                needOld = True
                # raise ValueError("Need AFlux and DFlux")

        elif isinstance(f, BndFlux_v):
            boundary_flux_vs.update([(kk, f) for kk in ids])

        elif isinstance(f, BndFlux_c):
            boundary_flux_cs.update([(kk, f) for kk in ids])

        elif isinstance(f, BndValue):
            boundary_values.update([(kk, f) for kk in ids])

        else:
            needOld = True
            # raise NotImplementedError(f"unknown boundary type {k} : {f}")

        #########################################################
        if needOld: # this will be removed
            deprecated("""Use new tag classes in Model's boundary dictionary:
    dune.femdg.BndValue|BndFlux_c|BndFlux_v""") # ,stacklevel=-1) # stacklevel=10 might work
            if isinstance(f,tuple) or isinstance(f,list):
                boundary_flux_cs.update( [ (kk,BndFlux_c(f[0])) for kk in ids] )
                boundary_flux_vs.update( [ (kk,BndFlux_v(f[1])) for kk in ids] )
            else:
                if len(inspect.signature(f).parameters) == 4:
                    if hasAdvFlux and not hasDiffFlux:
                        boundary_flux_cs.update( [ (kk,BndFlux_v(f)) for kk in ids] )
                    elif not hasAdvFlux and hasDiffFlux:
                        boundary_flux_vs.update( [ (kk,BndFlux_c(f)) for kk in ids] )
                    else:
                        assert not (hasAdvFlux and hasDiffFlux), "one boundary flux provided for id "+str(k)+" but two bulk fluxes given"
                else:
                    boundary_values.update( [ (kk,BndValue(f)) for kk in ids] )

    return boundary_flux_cs, boundary_flux_vs, boundary_values


def _splitBoundary(Model, override_boundary_dict=None):
    """ take a Model and split boundary dictionary testing and requirements are met
    """
    if override_boundary_dict is not None:
        boundary_dict = override_boundary_dict
    else:
        try:
            boundary_dict = Model.boundary
        except AttributeError:
            boundary_dict = {}

    hasFlux_c = hasattr(Model, "F_c")
    hasFlux_v = hasattr(Model, "F_v")

    boundary_flux_cs, boundary_flux_vs, boundary_values = classify_boundary(
        boundary_dict, hasFlux_c, hasFlux_v
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

def splitBoundary(Model, t,x,u,n):
    """ take a Model and split boundary dictionary testing and requirements are met
        Replaces all callables with the correct UFL expressions using
        t, SC, TrialF, FacetN
    """
    boundary_flux_cs, boundary_flux_vs, boundary_values = _splitBoundary(Model)

    boundary_flux_cs = {
        (k(x) if callable(k) else k): f(t, x, u, n) for k, f in boundary_flux_cs.items()
    }
    boundary_flux_vs = {
        (k(x) if callable(k) else k): f(t, x, u, ufl.grad(u), n)
        for k, f in boundary_flux_vs.items()
    }
    boundary_values = {
        (k(x) if callable(k) else k): f(t, x, u, n)
         for k, f in boundary_values.items()
    }
    hasBoundaryValue = {k: True for k in boundary_values.keys()}

    return (
        boundary_flux_cs,
        boundary_flux_vs,
        boundary_values,
        hasBoundaryValue,
    )
