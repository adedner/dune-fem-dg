import time, math, sys
from ufl import *
from dune.ufl import DirichletBC, Constant
from dune.grid import structuredGrid, cartesianDomain
from dune.alugrid import aluCubeGrid as cubeGrid
from dune.alugrid import aluSimplexGrid as simplexGrid
from dune.fem.space import lagrange, dgonb, raviartThomas
from dune.fem.scheme import galerkin
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper
import dune.create as create
from dune.grid import reader


def problem(V):

    class Model:
        # Set up time
        secperhour = 3600
        t_start_hour = 0
        t_start = t_start_hour * secperhour  # start (s)
        t_end_hour = 250  # (h)
        t_end = t_end_hour * secperhour  # (s)
        nt_hour = 36  # 12 time steps per hour
        nt = nt_hour * (t_end_hour - t_start_hour)  # time steps total
        dt_ = 10*(t_end - t_start) / nt  # time step (s)
        # t_values = np.arange(t_start, t_end + dt_, dt_)
        #t_values = np.linspace(t_start, t_end, nt)
        t = Constant(t_start, "time")  # time variable

        # reaction, diffusion and other coefficients
        # K   = dune.ufl.Constant(10.0, "reactionRate")
        # dt = dune.ufl.Constant(0.01, "timeStep")
        # t  = dune.ufl.Constant(0., "timeStep")
        phi_crit = 0.1
        # K = Constant(1.0e-12) * ((phi - phi_crit) / (phi0 - phi_crit)) ** 3  #
        rhoc = 2710  # density calcite
        Mb = 1.0  # Molar mass suspended biomass
        Mc = 0.1  # Molar mass calcite

        Mc_rhoc = Mc / rhoc

        #outflowVeldirich = conditional(x[0]>-0.3,1.,0.) * conditional(x[1]>-0.3,1.,0.) * conditional(x[0]<0.3,1.,0.) * conditional(x[1]<0.3,1.,0.)

        # Well
        Qp = 1.25e-3
        Qw = Qp
        # source_p = Well(mesh=mesh, Q=Qp, degree=1)  # UNCOMMENT
        # Pressure (natural)
        # p0_in = Constant(1.25e7)

        cu_in = 2000
        Qcu = Qp * cu_in

        cb_in = 0.01
        Qcb = Qp * cb_in
        #
        # Dispersivity (constant):
        #
        D_mol = 1.587e-9  # molecular diffusion
        alpha_long = 0.01  # longitudinal dispersivity
        # CONSTANT
        # D_mech = Constant(1e-6)
        D = D_mol

        # FULL TENSOR
        # I = Identity(dim)

        kub = 0.01
        kurease = 706.7
        ke = kub * kurease
        Ku = 355  # Ku = KU * rho_w = 0.355 * 1e3

        # attachment
        ka = 8.51e-7

        # fu=0
        # Biomass
        # From Johannes
        b0 = 3.18e-7  # decay rate coeff.
        c_attach = 4e-6  # unspecific attachment rate


        # Biomass
        # From Johannes
        b0 = 3.18e-7  # decay rate coeff.

        ## ufl expression for

        # pressure = [ -sin(2*pi*(x[0]-0.5)*(x[1]-0.5)) ]
        K = 1.0e-12

        # Misc. other parameters
        mu = 0.0008  # viscosity of water

        transportVelocity = V

        dimRange = 3 # three unknowns: (phi, phi cu, phi cb)^T

        dim = V.space.dimDomain # TODO extract from domain

        ### initial ###
        phi0 = 0.2
        cu0  = 0
        cb0  = 0

        ### Model functions ###
        def toPrim(U):
            # ( phi, phi cu, phi cb ) --> (phi, cu, cb)
            return U[0], U[1] / U[0], U[2] / U[0]

        def toCons(U):
            # ( phi, cu, cb ) --> (phi, phi cu, phi cb)
            return U[0], U[1] * U[0], U[2] * U[0]

        def darcyVelocity( p ):
            return -1./Model.mu * Model.K * grad(p[0])

        # right hand side for pressure equation
        def pressureRHS(U):
            phi, cu, cb = Model.toPrim( U )
            qu   = Model.ke * phi * Model.Mb * cb * cu / (Model.Ku + cu )
            hour = Model.t / Model.secperhour
            # TODO fix hours here
            Qw1_t = conditional( hour > 20., conditional( hour < 25., Model.Qw, 0), 0 )
            Qw2_t = conditional( hour > 45., conditional( hour < 60., Model.Qw, 0), 0 )
            return as_vector([ qu + Qw1_t + Qw2_t ])

        def S_s(t,x,U,DU): # or S_s for a stiff source
            phi, cu, cb = Model.toPrim( U )

            qu = Model.ke * phi * Model.Mb * cb * cu / (Model.Ku + cu )
            qc = qu # Assumption in the report, before eq 3.12
            qb = cb * ( phi * ( Model.b0 + Model.ka) + qc * Model.Mc_rhoc )

            hour = t / Model.secperhour
            Qu_t = conditional( hour > 25., conditional( hour < 45., Model.Qcu, 0), 0 )
            Qb_t = conditional( hour < 20., Model.Qcb, 0 )
            return as_vector([ -qu * Model.Mc_rhoc,
                               -qu + Qu_t,
                               -qb + Qb_t
                             ])
        def F_c(t,x,U):
            phi, cu, cb = Model.toPrim(U)
            # first flux should be zero since porosity is simply an ODE
            return as_matrix([ [*(Model.velocity(t,x,u) *  0)],
                               [*(Model.velocity(t,x,U) * cu)],
                               [*(Model.velocity(t,x,U) * cb)]
                             ])
        def maxLambda(t,x,U,n):
            return abs(dot(Model.velocity(t,x,U),n))
        def velocity(t,x,U):
            return Model.transportVelocity
        def F_v(t,x,U,DU):
            # eps * phi^(1/3) * grad U
            return Model.D * pow(U[0], 1./3.) * DU
        def maxDiffusion(t,x,U):
           return Model.D * pow(U[0], 1./3.)
        def physical(U):
            #phi, cu, cb = Model.toPrim( U )
            # U should be positive
            #return conditional( phi>=0,1,0)*conditional(cu>=0,1,0)*conditional(cb>=0,1,0)
            return conditional(U[0]>=0,1,0)*conditional(U[1]>=0,1,0)*conditional(U[2]>=0,1,0)
        def dirichletValue(t,x,u):
            return as_vector(Model.dimRange*[0])
        boundary = {(1,2,3,4): dirichletValue}

    Model.initial=as_vector([Model.phi0, Model.cu0, Model.cb0])
    Model.endTime=250 * Model.secperhour
    Model.name="micap"
    return Model
    # end def problem

#### main program ####

#gridView = cubeGrid(cartesianDomain([0,0],[1,1],[50,50]))
domain = (reader.gmsh, "circlemeshquad.msh")
gridView = create.view("adaptive", grid="ALUCube", constructor=domain, dimgrid=2)

# could also use dg here
pressureSpace = lagrange(gridView, order=1, dimRange=1)
pressure = pressureSpace.interpolate([1e-7],name="pressure")
rhs      = pressureSpace.interpolate([0],name="rhs")

u    = TrialFunction(pressureSpace)
phi  = TestFunction(pressureSpace)
x    = SpatialCoordinate(pressureSpace)
n    = FacetNormal(pressureSpace)


dbc  = DirichletBC(pressureSpace,[ 0 ])
# project into rt space
velocitySpace = raviartThomas(gridView,order=1,dimRange=1)

## projection seems to be buggy
# there interpolation is used
transportVelocity = velocitySpace.interpolate([0,0],name="velocity")

Model = problem( transportVelocity )

# create discrete function space
space = dgonb( gridView, order=2, dimRange=Model.dimRange)

# create and initialize solution
u_h = space.interpolate(Model.initial, name='u_h')

# TODO use Model.pressureRHS directly
#        Model.pressureRHS(u_h), phi) * dx
form = inner(grad(u),grad(phi)) * dx - inner(rhs, phi) * dx
pressureScheme = galerkin([form == 0, dbc], solver="cg",
                           parameters={"newton.linear.verbose":True,
                                       "newton.verbose":True} )

Stepper = femdgStepper(order=3, rkType="EX")

operator = femDGOperator(Model, space, limiter=None, threading=True)
operator.applyLimiter( u_h )
stepper  = Stepper(operator)

vtk = gridView.sequencedVTK("Chemical_highK", subsampling=1, celldata=[transportVelocity], pointdata=[pressure,u_h])
vtk() # output initial solution

# TODO t and Model.t
def updateVelocity( t ):
    Model.t.value = t
    # update right hand side
    rhs.interpolate( Model.pressureRHS( u_h ) )
    # re-compute pressure
    pressureScheme.solve(target=pressure)

    # Darcy velocity
    velocity = -1./Model.mu * Model.K * grad(pressure[0])

    # project into rt space
    transportVelocity.interpolate( velocity )

t        = 0
tcount   = 0
saveStep = 0.001 # Model.endTime/100
saveTime = saveStep
while t < Model.endTime:

    # compute new pressure
    # TODO, update pressure here

    print("### Compute Darcy velocity ###")
    updateVelocity( t )

    vtk()

    print("### Compute advection-diffusion ###")
    operator.setTime(t)
    dt = stepper(u_h)
    t += dt
    tcount += 1
    # check that solution is meaningful
    if math.isnan( u_h.scalarProductDofs( u_h ) ):
        vtk()
        print('ERROR: dofs invalid t =', t,flush=True)
        print('[',tcount,']','dt = ', dt, 'time = ',t, flush=True )
        sys.exit(1)
    if 1: # t > saveTime:
        # TODO: issue with time step estimate: here the 'fullPass' is used
        # but that is not set for IMEX
        print('[',tcount,']','dt = ', dt, 'time = ',t,
                'dtEst = ',operator.timeStepEstimate,
                'elements = ',Model.domain.size(0), flush=True )
        vtk()
        saveTime += saveStep

vtk() # output final solution
