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

# eps = 0.2                 # diffusion rate
K = 200                    # reaction rate
Q = 1                     # source strength
P1 = as_vector([0.1,0.1]) # midpoint of first source
P2 = as_vector([0.9,0.9]) # midpoint of second source
RF = 0.075                # radius of sources

def problem():

    def computeVelocityCurl():
        gridView = cubeGrid(cartesianDomain([0,0],[1,1],[20,20]))
        # TODO without dimRange the 'vectorization' is not carried out correctly
        streamSpace = lagrange(gridView, order=1, dimRange=1)
        Psi  = streanSpace.interpolate(0,name="streamFunction")
        u    = TrialFunction(streamSpace)
        phi  = TestFunction(streamSpace)
        x    = SpatialCoordinate(streamSpace)
        form = ( inner(grad(u),grad(phi)) -
                 0.1*2*(2*pi)**2*sin(2*pi*x[0])*sin(2*pi*x[1]) * phi[0] ) * dx
        dbc  = DirichletBC(streamSpace,[0])
        streamScheme = galerkin([form == 0, dbc], solver="cg")
        streamScheme.solve(target=Psi)
        velocity = as_vector([-Psi[0].dx(1),Psi[0].dx(0)]) # note: not extra projection
        return gridView, velocity

    def computeVelocityGrad():
        gridView = cubeGrid(cartesianDomain([0,0],[1,1],[50,50]))
        # could also use dg here
        pressureSpace = lagrange(gridView, order=1, dimRange=1)
        pressure      = pressureSpace.interpolate(0,name="pressure")
        u    = TrialFunction(pressureSpace)
        phi  = TestFunction(pressureSpace)
        x    = SpatialCoordinate(pressureSpace)
        n    = FacetNormal(pressureSpace)
        form = inner(grad(u),grad(phi)) * dx
        dbc  = DirichletBC(pressureSpace,[ -sin(2*pi*(x[0]-0.5)*(x[1]-0.5)) ])
        pressureScheme = galerkin([form == 0, dbc], solver="cg",
                                   parameters={"newton.linear.verbose":True,
                                               "newton.verbose":True} )
        pressureScheme.solve(target=pressure)
        ## ufl expression for

        # pressure = [ -sin(2*pi*(x[0]-0.5)*(x[1]-0.5)) ]
        K = Constant(1.0e-12)

        # Misc. other parameters
        mu = Constant(0.0008)  # viscosity of water

        # Darcy velocity
        velocity = -1./mu * K * grad(pressure[0])

        # project into rt space
        velocitySpace = raviartThomas(gridView,order=1,dimRange=1)
        ## projection seems to be buggy
        _velo = velocitySpace.project(velocity,name="velocity")
        velo = velocitySpace.interpolate(velocity,name="velocity")
        gridView.writeVTK("velocity",subsampling=1,
                pointvector={"pVelo":velo, "pProjVelo":_velo},
                cellvector={"cVelo":velo, "cProjVelo":_velo}
                )

        return gridView, velo

    class Model:
        # domain, transportVelocity = computeVelocityCurl()

        domain, transportVelocity = computeVelocityGrad()
        dimRange = 3 # three unknowns: (phi, phi cu, phi cb)^T

        dim = 2 # TODO extract from domain

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
        t = Constant(t_start)  # time variable

        # reaction, diffusion and other coefficients
        # K   = dune.ufl.Constant(10.0, "reactionRate")
        eps = Constant(0.01, "diffusionRate")
        # dt = dune.ufl.Constant(0.01, "timeStep")
        # t  = dune.ufl.Constant(0., "timeStep")
        phi_crit = Constant(0.1)
        # K = Constant(1.0e-12) * ((phi - phi_crit) / (phi0 - phi_crit)) ** 3  #
        rhoc = Constant(2710)  # density calcite
        Mb = Constant(1.0)  # Molar mass suspended biomass
        Mc = Constant(0.100)  # Molar mass calcite

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
        D_mol = Constant(1.587e-9)  # molecular diffusion
        alpha_long = 0.01  # longitudinal dispersivity
        # CONSTANT
        D_mech = Constant(1e-6)
        D = Constant(1e-6)

        # FULL TENSOR
        # I = Identity(dim)

        kub = 0.01
        kurease = 706.7
        ke = Constant(kub * kurease)
        Ku = Constant(355)  # Ku = KU * rho_w = 0.355 * 1e3

        # fu=0
        # Biomass
        # From Johannes
        b0 = Constant(3.18e-7)  # decay rate coeff.
        c_attach = Constant(4e-6)  # unspecific attachment rate


        # Biomass
        # From Johannes
        b0 = Constant(3.18e-7)  # decay rate coeff.

        def S_s(t,x,U,DU): # or S_s for a stiff source
            phi  = U[0] # porosity
            qu   = kurease * kub * phi * Mb * U[2] * U[1] / (Ku + U[1] )
            qc = qu # Assumption in the report, before eq 3.12
            qb = U[2] * ( phi * ( b0 + ka) + qc * Mc / rhoc )

            Qu_t = conditional( t > 25., conditional( t < 45., Qu, 0), 0 )
            Qb_t = conditional( t < 20., Qb, 0 )
            return as_vector([ -Mc * qu / rhoc,
                               -qu + Qu_t,
                               -qb + Qb_t
                             ])
        def F_c(t,x,U):
            # first flux should be zero since porosity is simply an ODE
            return as_matrix([  as_vector([ 0 ]),
                              *([Model.velocity(t,x,U) * U[i]/U[0] for i in range(1,dimRange)]),
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
            # U should be positive
            return conditional(U[0]>=0,1,0)*conditional(U[1]>=0,1,0)*conditional(U[2]>=0,1,0)
        def dirichletValue(t,x,u):
            return as_vector(Model.dimRange*[0])
        boundary = {(1,2,3,4): dirichletValue}

    Model.initial=as_vector([0,0,0])
    Model.endTime=250 * Model.secperhour
    Model.name="micap"
    return Model

Model = problem()
Stepper = femdgStepper(order=3)

# create discrete function space
space = dgonb( Model.domain, order=3, dimRange=Model.dimRange)
# operator = femDGOperator(Model, space, limiter="scaling", threading=True)
operator = femDGOperator(Model, space, limiter=None, threading=True)
stepper  = Stepper(operator)
# create and initialize solution
u_h = space.interpolate(Model.initial, name='u_h')
operator.applyLimiter( u_h )

vtk = Model.domain.sequencedVTK("Chemical_highK", subsampling=1, pointdata=[u_h])
vtk() # output initial solution

t        = 0
tcount   = 0
saveStep = 0.001 # Model.endTime/100
saveTime = saveStep
while t < Model.endTime:
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
