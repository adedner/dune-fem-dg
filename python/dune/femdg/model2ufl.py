from ufl import (
    FacetNormal,
    SpatialCoordinate,
    TestFunction,
    TrialFunction,
    as_vector,
    ds,
    dx,
    grad,
    inner,
)

try:
    import dune.ufl as duneUfl
except:
    duneUfl = None
from .boundary import splitBoundary

# a simple class used in case dune.ufl.DirichletBC could not be used, i.e.,
# some other discretization package is used
class DefaultDirichletBC:
    def __init__(self, space, value, domain=None):
        self.space = space
        self.value = value
        self.domain = domain

    def __str__(self):
        return str(self.value) + str(self.domain)

class DefaultConstant:
    def __init__(self, value, name):
        self.value = value
        self.name = name

def model_ufl(Model, space, t, DirichletBC, Constant):
    u = TrialFunction(space)
    v = TestFunction(space)
    x = SpatialCoordinate(space.cell())
    n = FacetNormal(space.cell())

    f_c_model = None
    if hasattr(Model, "F_c"):
        f_c_model = inner(Model.F_c(t, x, u), grad(v)) # -div F_c v
    if hasattr(Model, "S_e"):
        # there is an issue with S_e returning 'zero' and zero*dx leading to UFL error
        se = (
            inner(as_vector(Model.S_e(t, x, u, grad(u))), v)
        )  # (-div F_c + S_e) * v
        if f_c_model is not None:
            f_c_model += se
        else:
            f_c_model = se

    f_v_model = None
    if hasattr(Model, "F_v"):
        f_v_model = inner(Model.F_v(t, x, u, grad(u)), grad(v)) # -div F_v v

    if hasattr(Model, "S_i"):
        si = inner(as_vector(Model.S_i(t, x, u, grad(u))), v) # (-div F_v + S_i) v
        if f_v_model is not None:
            f_v_model += si
        else:
            f_v_model = si

    # need to extract boundary information from Model
    (
        boundary_flux_cs,
        boundary_flux_vs,
        boundary_values,
        hasBoundaryValue,
    ) = splitBoundary(Model, t, x, u, n)

    dirichletBCs = [
        DirichletBC(space, item[1], item[0]) for item in boundary_values.items()
    ]
    boundary_flux_vs = -sum(
        [inner(item[1], v) * ds(item[0]) for item in boundary_flux_vs.items()]
    )  # keep all forms on left hand side
    boundary_flux_cs = -sum(
        [inner(item[1], v) * ds(item[0]) for item in boundary_flux_cs.items()]
    )  # keep all forms on left hand side

    if f_c_model:
        f_c_model = f_c_model * dx
    if f_v_model:
        f_v_model = f_v_model * dx
    # !!! fix issue with f_?_model==zero not being a form
    return (
        f_c_model,
        f_v_model,
        {
            "dirichletBCs": dirichletBCs,
            "boundary_flux_cs": boundary_flux_cs,
            "boundary_flux_vs": boundary_flux_vs,
            "hasBoundaryValue": hasBoundaryValue,
        },
    )


def model2ufl(
    Model, space, initialTime=0, *,
    DirichletBC=DefaultDirichletBC, Constant=DefaultConstant,
    returnFull=False
):
    if duneUfl is not None:
        try:
            dimDomain = space[0]
            dimRange = space[1]
            space = duneUfl.Space(dimDomain,dimRange=dimRange)
        except:
            pass
        if DirichletBC == DefaultDirichletBC:
            DirichletBC = duneUfl.DirichletBC
        if Constant == DefaultConstant:
            Constant = duneUfl.Constant
    class M(Model):
        if hasattr(Model, "S_e"):

            def S_e(t, x, U, DU):
                return -Model.S_e(t, x, U, DU)

        if hasattr(Model, "S_i"):

            def S_i(t, x, U, DU):
                return -Model.S_i(t, x, U, DU)

        if hasattr(Model, "F_c"):

            def F_c(t, x, U):
                return -Model.F_c(t, x, U)

    timeC = Constant(initialTime,"time")
    f_c_model, f_v_model, boundary_model = model_ufl(M, space, timeC, DirichletBC, Constant)
    boundary_model["boundary_flux_cs"] = -boundary_model["boundary_flux_cs"]
    form = boundary_model["boundary_flux_cs"] + boundary_model["boundary_flux_vs"]
    if f_c_model is not None:
        form += f_c_model
    if f_v_model is not None:
        form += f_v_model

    if not returnFull:
        return [form == 0, *boundary_model["dirichletBCs"]]
    else:
        boundary_model["f_c_model"] = f_c_model
        boundary_model["f_v_model"] = f_v_model
        boundary_model["form"] = form
        return boundary_model

# full set of boundary conditions still need to be added
def model2dgufl(Model,space):
    from ufl import ( TrialFunction, TestFunction,
                      SpatialCoordinate, FacetNormal,
                      dx, ds, dot )
    from dune.femdg.boundary import splitBoundary
    from dune.ufl import Constant
    from dune.femdg.dolfin_dg import (
            HyperbolicOperator, DGDirichletBC, LocalLaxFriedrichs,
            EllipticOperator, DGFemNIPG, DGFemBO, DGFemSIPG
        )
    t = Constant(0,"time")
    x = SpatialCoordinate(space.cell())
    n = FacetNormal(space.cell())
    u = TrialFunction(space)
    v = TestFunction(space)

    ( boundary_flux_cs,
      boundary_flux_vs,
      boundary_values,
      hasBoundaryValue ) = splitBoundary(Model, t, x, u, n)

    # is there a way to get different boundary condition into the dolfin-dg operators?
    # assume single dirichlet condition for now.
    # Perhaps ds(id) can work?
    dbc = []
    for i,item in enumerate(boundary_values.items()):
        dbc += [DGDirichletBC(ds(item[0]), item[1])]

    rhs = 0
    if hasattr(Model,"S_i"):
        rhs += dot(Model.S_i(t,x,u,grad(u)),v) * dx
    if hasattr(Model,"S_e"):
        rhs += dot(Model.S_e(t,x,u,grad(u)),v) * dx

    lhs = 0
    if hasattr(Model,"F_c"):
        def F_c(u):
            return Model.F_c(t,x,u)
        if hasattr(Model,"maxWaveSpeed"):
            def alpha(u, n):
                return Model.maxWaveSpeed(t,x,u,n)
        else:
            def alpha(u, n):
                return Model.dimRange*[0.]
        ho = HyperbolicOperator(space.cell(), space, dbc, F_c,
                                LocalLaxFriedrichs(alpha))
        lhs += ho.generate_fem_formulation(u, v)
    if hasattr(Model,"F_v"):
        def F_v(u, grad_u):
            return Model.F_v(t,x,u,grad_u)
        eo = EllipticOperator(space.cell(), space, dbc, F_v)
        lhs += eo.generate_fem_formulation(u, v, vt = DGFemSIPG)

    return lhs == rhs
