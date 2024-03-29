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

A detailed description of how to use the code can be found in the first paper,
the description of the schemes can be found
in the second and third papers and an overview on performance of
the code is given in the forth paper.

* A. Dedner, R. Klöfkorn.
Extendible and Efficient Python Framework for Solving Evolution Equations
with Stabilized Discontinuous Galerkin Method. Commun. Appl. Math. Comput., 2021. https://dx.doi.org/10.1007/s42967-021-00134-5.

* S. Brdar, A. Dedner, and R. Klöfkorn.
Compact and stable Discontinuous Galerkin methods for convection-diffusion problems.
SIAM J. Sci. Comput., 34(1):263-282, 2012. https://dx.doi.org/10.1137/100817528

* A. Dedner and R. Klöfkorn. A Generic Stabilization Approach for Higher Order Discontinuous Galerkin Methods for Convection Dominated Problems.
J. Sci. Comput., 47(3):365-388, 2011. https://dx.doi.org/10.1007/s10915-010-9448-0

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

**Using pip**

dune-fem-dg can be installed using the Package Index of Python (pip).

```
pip install dune-fem-dg
```

See https://dune-project.org/doc/installation-pip/ for a more detailed
description.

**From source**

dune-fem-dg also provides a [shell script](https://gitlab.dune-project.org/dune-fem/dune-fem-dg/-/blob/master/scripts/build-dune-fem-dg.sh)
to install the sources using a specific git branch.
Further detailed explanation on the DUNE installation process please read
the installation notes https://www.dune-project.org/doc/installation/.

Documentation
-------------

A documentation of the Python bindings of DUNE-FEM-DG can be found in

* A. Dedner, R. Klöfkorn.
Extendible and Efficient Python Framework for Solving Evolution Equations
with Stabilized Discontinuous Galerkin Method. Commun. Appl. Math. Comput., 2021. https://dx.doi.org/10.1007/s42967-021-00134-5.

A documentation for C++ code (and 2.4 release) can be found in

* A. Dedner, S. Girke, R. Klöfkorn, T. Malkmus. The DUNE-FEM-DG module.
Archive of Numerical Software 5(1): 21--61. 2017. https://dx.doi.org/10.11588/ans.2017.1.28602
