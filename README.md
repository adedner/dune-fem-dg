DUNE-FEM-DG
===========

DUNE-FEM-DG is the implementation of Discontinuous Galerkin schemes
using the DUNE-FEM framework. Stabilized DG schemes for hyperbolic
as well as a wide range of different primal formulations
for elliptic/parabolic problems are implemented.
The operators can be used efficiently both in a
explicit/matrix free implementation or used to
setup a system matrix for use with the linear solvers available in DUNE-FEM.


License
-------

The DUNE-FEM-DG module is available under
the GNU General Public License version 2, or (at your option),
any later version.


References
----------

A detailed description of the schemes can be found
in the first two papers and an overview on performance of
the code is given in the third paper.

* S. Brdar, A. Dedner, and R. Klöfkorn.
Compact and stable Discontinuous Galerkin methods for convection-diffusion problems.
SIAM J. Sci. Comput., 34(1):263-282, 2012. http://dx.doi.org/10.1137/100817528

* A. Dedner and R. Klöfkorn. A Generic Stabilization Approach for Higher Order Discontinuous Galerkin Methods for Convection Dominated Problems.
J. Sci. Comput., 47(3):365-388, 2011. http://dx.doi.org/10.1007/s10915-010-9448-0

An overview on performance of the code is given in

* R. Klöfkorn. Efficient Matrix-Free Implementation of Discontinuous Galerkin Methods for Compressible Flow Problems.
Proceedings of the ALGORITMY 2012. http://www.iam.fmph.uniba.sk/algoritmy2012/zbornik/2Kloefkornf.pdf

By using the code you agree to cite one or both of the first two papers in any publication using this code.


Eye-candy
---------

The avatar of the project shows the solution of
the compressible Euler equations in 3D using
the parallel-adaptive DUNE-ALUGrid and the DG
discretization implemented in DUNE-FEM-DG.


Installation
------------

For a full explanation of the DUNE installation process please read
the installation notes https://www.dune-project.org/doc/installation/.


Documentation
-------------

A documentation for the 2.4 release can be found in

* A. Dedner, S. Girke, R. Klöfkorn, T. Malkmus. 2017. The DUNE-FEM-DG module. 
Archive of Numerical Software 5(1): 21--61. http://dx.doi.org/10.11588/ans.2017.1.28602
