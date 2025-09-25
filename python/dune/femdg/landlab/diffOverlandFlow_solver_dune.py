#!/usr/bin/env python3
"""Solve advection numerically using dune-fem-dg module."""

import numpy as np
from matplotlib import pyplot as plt

# import dolfin_dg as dg
import dune.femdg.dolfin_dg as dg

from landlab import Component
from landlab import LinkStatus
from landlab.field.errors import FieldError
from landlab.utils.return_array import return_array_at_link
from landlab.utils.return_array import return_array_at_node

from ufl import as_vector, as_matrix, exp, grad, div, conditional, dot, zero, sqrt
from dune.ufl import cell, Constant

from dune.fem.space import dgonb, dglegendre
from dune.fem.scheme import galerkin
from dune.fem.plotting import plotPointData as plot
from dune.femdg import model2ufl, model2dgufl
from dune.femdg import BndValue, BndFlux_v, BndFlux_c
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper

from .dune_landlab_adapter import DuneLandlabAdapter

class DuneDiffOverlandFlowSolver(Component):
    ## class variables

    _name = "DuneDiffOverlandFlowSolver"

    _unit_agnostic = True

    _info = {
        "surface_water__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of water on the surface",
        },
        "surface_water__depth_at_link": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "link",
            "doc": "Depth of water on the surface at grid links",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water__specific_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m2/s",
            "mapping": "link",
            "doc": "flow discharge component in the direction of the link",
        },
        "water__velocity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "flow velocity component in the direction of the link",
        },
        "water_surface__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "Downstream gradient of the water surface.",
        },
        "water_surface__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Elevation of the water surface.",
        },
    }

    def model(self):
        class ModelH:
            dimRange = 1
            def Hpos(U):
                if self._linear:
                    # return self.uhOld[0]
                    return conditional(self.uhOld[0]>0,self.uhOld[0],0)
                else:
                    # return U[0]
                    return conditional(U[0]>0,U[0],0)
            def D(U):
                Hpos = ModelH.Hpos(U)
                return Hpos**(7/3) / (self._n**2*self._uc)
            def F_v(t,x,U,DU):
                d = ModelH.D(U)
                dw = DU[0,:] + self.dEta
                f = d * as_vector([[dw[0]*self.xiRT[0],dw[1]*self.xiRT[1]]])
                f += self._eps*DU
                return f
            def S_i(t,x,U,DU):
                xi = self.xi0[0]
                Hpos = ModelH.Hpos(U)
                I = self._Ic * (1-exp(-Hpos/self._Hi))
                S = xi * as_vector([self._R - I])
                S -= (U - self.uhOld) / self._dt
                return S
            boundary = {range(1,5): BndValue(lambda t,x,U: as_vector([0]))}
        return ModelH

    def __init__(
        self,
        grid,
        *,
        roughness=0.01,               # n
        rain_rate=1.0e-5,             # R
        infilt_rate=0.0,              # Ic
        infilt_depth_scale=0.001,     # Hi
        velocity_scale=1.0,           # uc
        method="linear",
        eps=0.1
    ):
        self._linear = (method == "linear")
        # Call base class methods to check existence of input fields,
        # create output fields, etc.
        super().__init__(grid)
        self.initialize_output_fields()
        self._depth = grid.at_node["surface_water__depth"]
        self._depth_at_link = grid.at_link["surface_water__depth_at_link"]
        self._vel = grid.at_link["water__velocity"]
        self._disch = grid.at_link["water__specific_discharge"]
        self._wsgrad = grid.at_link["water_surface__gradient"]
        self._water_surf_elev = grid.at_node["water_surface__elevation"]
        self._inactive_links = grid.status_at_link == grid.BC_LINK_IS_INACTIVE

        self._R  = Constant(rain_rate,"rainRate")
        self._Ic = Constant(infilt_rate,"maxInfiltationRate")
        self._Hi = Constant(infilt_depth_scale,"charactWaterDepth")
        self._n  = Constant(roughness,"Manning")
        self._uc = Constant(velocity_scale,"charactScaleVelo")
        self._dt = Constant(1.,"dt")
        self._eps = Constant(eps,"eps")

        self._landlabAdapter = DuneLandlabAdapter(grid)
        self.gridView = self._landlabAdapter.gridView

        # get topography and its gradients
        self._elev = grid.at_node["topographic__elevation"]
        gradElev = grid.calc_grad_at_link(self._elev)
        gradElev[grid.status_at_link == grid.BC_LINK_IS_INACTIVE] = 0.0

        # setup topography as FV and Lagrange function

        self.eta0 = self._landlabAdapter.nodeFct(name="elev0")
        self._landlabAdapter.fromNode(self._elev, self.eta0)
        self.eta1 = self._landlabAdapter.cellFct(name="elev1")
        self.eta1.project(self.eta0)
        # setup gradient of eta as RT function
        self.dEta = self._landlabAdapter.linkFct(name="gradElev")
        self._landlabAdapter.fromLink(gradElev,self.dEta)

        # setup characteristic function for the basin (node and link)
        self.xi0 = self._landlabAdapter.nodeFct(name="basinNode")
        self._landlabAdapter.fromNode(self._elev, self.xi0)
        self.xi0.as_numpy[np.isclose(self._elev, -9999.0)] = 0
        self.xi0.as_numpy[self._elev>-9000] = 1
        gradElev[:] = 0
        gradElev[grid.status_at_link == grid.BC_LINK_IS_ACTIVE] = 1
        self.xiRT = self._landlabAdapter.linkFct(name="basinLink")
        self._landlabAdapter.fromLink(gradElev,self.xiRT)

        self.gridView.writeVTK("elev",celldata=[self.xi0,self.xiRT,self.dEta])

        self.Model = self.model()
        self.uh = dgonb(self.gridView,dimRange=self.Model.dimRange,order=1).function(name="uh")
        # self.uh = dglegendre(self.gridView,dimRange=1,order=1).function(name="uh")
        self.uhOld = self.uh.copy(name="old")
        self._advector = galerkin( model2dgufl(self.Model, self.uh.space),
                   solver="gmres", # "cg",
                   parameters={"linear.preconditioning.method":"ssor",   # none < ssor < ilu...
                               "linear.verbose":False,
                               "nonlinear.verbose":False,
                   })
        self._basisSize = self.uh.space.localBlockSize

    @property
    def rain_rate(self):
        return self._R.value
    @rain_rate.setter
    def rain_rate(self,value):
        self._R.value = value

    def _toLandlab(self, dh):
        self._depth[:] = dh
        self._depth.clip(min=0.0, out=self._depth)
        self._water_surf_elev[:] = self._elev + self._depth
        self.grid.map_value_at_max_node_to_link(
                        self._water_surf_elev, "surface_water__depth",
                        out=self._depth_at_link
                   )
        self.grid.calc_grad_at_link(self._water_surf_elev, out=self._wsgrad)
        self._wsgrad[self._inactive_links] = 0.0
        self._vel[:] = ( - self._depth_at_link**(4/3) * self._wsgrad /
                                             (self._uc.value*self._n.value**2) )
        self._disch[:] = self._depth_at_link * self._vel

    def update(self, dt):
        """Update the solution by one time step dt.

        Same as :meth:`~.run_one_step`.

        Parameters
        ----------
        dt : float
            Time-step duration. Needed to calculate the Courant number.
        """

        self._dt.value = dt
        self.uhOld.assign(self.uh)
        info = self._advector.solve(target=self.uh)
        self._toLandlab(self.uh.as_numpy[0::self._basisSize])

    def run_one_step(self, dt):
        """Update the solution by one time step dt.

        Same as :meth:`~.update`.

        Parameters
        ----------
        dt : float
            Time-step duration. Needed to calculate the Courant number.
        """
        self.update(dt)

class DuneDiffOverlandFlowSolver1(DuneDiffOverlandFlowSolver):
    _name = "DuneDiffOverlandFlowSolver1"
    def __init__(
        self,
        grid,
        *,
        roughness=0.01,               # n
        rain_rate=1.0e-5,             # R
        infilt_rate=0.0,              # Ic
        infilt_depth_scale=0.001,     # Hi
        velocity_scale=1.0,           # uc
        method="linear",
        eps=0.1
    ):
        super().__init__(grid,roughness=roughness,rain_rate=rain_rate,infilt_rate=infilt_rate,
                              infilt_depth_scale=infilt_depth_scale,
                              velocity_scale=velocity_scale,method=method,eps=eps)
        self.uh.interpolate([0,self.dEta[0],self.dEta[1]])
    def model(self):
        class Model:
            # original  S(H) - H_t = - div( H^(7/3) / (self._n**2*self._uc) grad(H + eta) )
            # modified  S(H) - H_t = - div( H^(4/3) / (self._n**2*self._uc) W H)
            #                    W = grad(H + eta)
            dimRange = 3
            def Hpos(U):
                return conditional(U[0]>0,U[0],0)
            def velocity(U):
                H = Model.Hpos(U)
                c = 1 / (self._n**2*self._uc)
                return c * as_vector( [U[1]*H**(4/3), U[2]*H**(4/3)] )
            def F_c(t,x,U):
                H = U[0]
                vel = Model.velocity(U)
                f = as_matrix([ [ -self.xiRT[0] * vel[0]*H , -self.xiRT[1] * vel[1]*H ],
                                [ -self.xiRT[0] * H        , 0                        ],
                                [ 0                        , -self.xiRT[1] * H        ] ])
                return f
            def maxWaveSpeed(t,x,U,n):
                return U[1]*n[0] + U[2]*n[1]
            def S_i(t,x,U,DU):
                H = Model.Hpos(U)
                I = self._Ic * (1-exp(-H/self._Hi))
                S = as_vector([ self.xi0[0] * (self._R - I), 0, 0 ])
                S -= (U - self.uhOld) / self._dt
                return S
            boundary = {range(1,5): BndValue(lambda t,x,U: as_vector([0,0,0]))}
        return Model

    def update(self, dt):
        """Update the solution by one time step dt.

        Same as :meth:`~.run_one_step`.

        Parameters
        ----------
        dt : float
            Time-step duration. Needed to calculate the Courant number.
        """

        self._dt.value = dt
        self.uhOld.assign(self.uh)
        info = self._advector.solve(target=self.uh)
        print(info)
        self._toLandlab(self.uh.as_numpy[0::self._basisSize])
