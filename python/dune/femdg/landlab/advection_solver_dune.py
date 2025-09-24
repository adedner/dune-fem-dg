#!/usr/bin/env python3
"""Solve advection numerically using dune-fem-dg module."""

import numpy as np

from landlab import Component
from landlab import LinkStatus
from landlab.field.errors import FieldError
from landlab.utils.return_array import return_array_at_link
from landlab.utils.return_array import return_array_at_node

### dune imports ###
from dune.grid import cartesianDomain
from dune.grid import yaspGrid as gridView

from dune.fem.space import dgonb, finiteVolume, raviartThomas, dglegendre
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper
from ufl import SpatialCoordinate, as_matrix, as_vector, dot, conditional
from dune.ufl import cell
from dune.femdg import BndValue, BndFlux_v, BndFlux_c

from .dune_landlab_adapter import DuneLandlabAdapter

class DuneAdvectionSolver(Component):
    """Numerical solution for advection using a higher order finite volume method based on dune-fem-dg.

    The component is restricted to regular grids (e.g., Raster or Hex).
    If multiple fields are advected, the advection__flux field will apply to
    the last one listed.

    Parameters
    ----------
    grid : RasterModelGrid or HexModelGrid
        A Landlab grid object.
    fields_to_advect : field name or list or (n_nodes,) array (default None)
        A node field of scalar values that will be advected, or list of fields.
        If not given, the component creates a generic field, initialized to zeros,
        called advected__quantity. If list >1 element given, advection will be
        applied to each field in it.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import DuneAdvectionSolver
    >>> grid = RasterModelGrid((3, 7))
    >>> s = grid.add_zeros("advected__quantity", at="node")
    >>> s[9:12] = np.array([1.0, 2.0, 1.0])
    >>> u = grid.add_zeros("advection__velocity", at="link")
    >>> u[grid.horizontal_links] = 1.0
    >>> advec = DuneAdvectionSolver(grid, fields_to_advect="advected__quantity")
    >>> for _ in range(5):
    ...     advec.update(0.2)
    ...
    >>> np.argmax(s[7:14])
    4
    """

    ## class variables

    _name = "DuneAdvectionSolver"
    _counter = 0

    _unit_agnostic = True

    _info = {
        "advected__quantity": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Scalar quantity advected",
        },

        "advection__velocity": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/y",
            "mapping": "link",
            "doc": "Link-parallel advection velocity magnitude",
        },
    }

    def advectionModel(self,dim, dimRange, velocity):
        x = SpatialCoordinate(cell(dim))
        eps = None # no diffusion
        dimR = dimRange

        v = velocity

        class Model:
            dimRange = dimR ## number of fields to advect
            if v is not None:
                # advective flux
                def F_c(t,x,U):
                    return as_matrix([ [*(v*u)] for u in U ])

                # max wave speed to compute dt
                def maxWaveSpeed(t,x,U,n):
                    return abs(dot(v,n))
                # velocity used in advection
                def velocity(t,x,U):
                    return v

            if eps is not None:
                def F_v(t,x,U,DU):
                    return eps*DU
                def maxDiffusion(t,x,U):
                   return eps

            lowerBound = [0]*dimRange
            upperBound = [1]*dimRange

            # simple 'dirchlet' boundary conditions on all boundaries
            boundary = {range(1,5): BndValue( lambda t,x,U: as_vector([0]*dimRange) )}

        return Model

    def __init__(
        self,
        grid,
        fields_to_advect=None,
        *,
        limiter="default",  # default is 'scalling' with given model
        initial=None,
        order=0,            # 0: FV, 'order'>0: dg
        velocity=None,      # None: use velocity on grid links, otherwise assume 'velocity' is a list/tuple
        **kwargs
    ):
        """Initialize DuneAdvectionSolver."""

        # Call base class methods to check existence of input fields,
        # create output fields, etc.
        super().__init__(grid)
        self.initialize_output_fields()

        self._scalars = []  # list of fields to advect
        self._fluxes = []  # list of flux fields
        if fields_to_advect is None:
            try:
                self._scalars.append(self.grid.at_node["advected__quantity"])
            except KeyError:
                self._scalars.append(
                    self.grid.add_zeros("advected__quantity", at="node")
                )
            try:
                self._fluxes.append(self.grid.at_link["advection__flux"])
            except KeyError:
                self._fluxes.append(self.grid.add_zeros("advection__flux", at="link"))
        elif isinstance(fields_to_advect, list):
            flux_counter = 0
            for field in fields_to_advect:
                self._scalars.append(return_array_at_node(self.grid, field))
                if isinstance(field, str):
                    flux_name = "flux_of_" + field
                else:
                    flux_name = "advection__flux_" + str(flux_counter)
                    flux_counter += 1
                try:
                    flux = return_array_at_link(self.grid, flux_name)
                except FieldError:
                    flux = grid.add_zeros(flux_name, at="link")
                self._fluxes.append(flux)
        else:
            self._scalars.append(return_array_at_node(self.grid, fields_to_advect))
            if isinstance(fields_to_advect, str):
                flux_name = "flux_of_" + fields_to_advect
            else:
                flux_name = "advection__flux"
            try:
                flux = return_array_at_link(self.grid, flux_name)
            except FieldError:
                flux = grid.add_zeros(flux_name, at="link")
            self._fluxes.append(flux)

        self._vel = self.grid.at_link["advection__velocity"]

        # Setup dune structures
        self._landlabAdapter = DuneLandlabAdapter(grid)
        self._gridView = self._landlabAdapter.gridView
        dimRange = len(self._scalars)

        if order is None or order == 0:
            self._uh = self._landlabAdapter.fromNode(self._scalars, "uh")
            self.dg = False
        else:
            self._space = dgonb( self._gridView, order=order, dimRange=dimRange)
            self.dg = True
            self._uh = self._space.function(name="uh")

        if velocity is None:
            self._velocity = self._landlabAdapter.linkFct("velocity"+str(DuneAdvectionSolver._counter))
            self._landlabAdapter.fromLink(self._vel,self._velocity)
            DuneAdvectionSolver._counter += 1
        else:
            self._velocity = as_vector(velocity)

        # number of scalar fields to advect
        AdvModel = self.advectionModel( self._gridView.dimGrid, dimRange=dimRange,
                                        velocity=self._velocity
                                       )

        self.initial(initial)

        codegen    = False
        threading  = False

        ## create fem-dg operator
        self._op = femDGOperator(AdvModel, self._uh.space,
                                 limiter=limiter,
                                 advectionFlux="Dune::Fem::AdvectionFlux::Enum::upwind",
                                 threading=threading,
                                 codegen=codegen)

        # time stepping order of RK scheme
        rkOrder = self._uh.space.order + 1

        ## get class
        Stepper = femdgStepper(order=rkOrder)

        ## create stepper object
        self._advector = Stepper( self._op, cfl=1.0 )

    def initial(self,initial=None):
        if initial:
            self._uh.interpolate(initial)
            self.reinitial = False
        else:
            self._landlabAdapter.fromNode( self._scalars, self._uh )
            # set from outside to allow changing the averages from outside
            # but that then we need to call limiter by hand as well
            self.reinitial = not self.dg

    def update(self, dt):
        """Update the solution by one time step dt.

        Same as :meth:`~.run_one_step`.

        Parameters
        ----------
        dt : float
            Time-step duration. Needed to calculate the Courant number.
        """

        ## update velocity field
        if hasattr(self._velocity,"space"):
            self._landlabAdapter.fromLink(self._vel, self._velocity)

        dimRange = len(self._scalars)
        ## self._scalars --> self._uh
        if self.reinitial:
            self._landlabAdapter.fromNode( self._scalars, self._uh )

        # advect by fixed dt
        self._advector(self._uh, dt )

        ## copy back to each component (or average)
        self._landlabAdapter.toNode( self._uh, self._scalars )

    def run_one_step(self, dt):
        """Update the solution by one time step dt.

        Same as :meth:`~.update`.

        Parameters
        ----------
        dt : float
            Time-step duration. Needed to calculate the Courant number.
        """
        self.update(dt)
