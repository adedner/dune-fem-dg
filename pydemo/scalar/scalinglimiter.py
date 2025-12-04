import matplotlib.pyplot as plt
import numpy as np

try:
    from dune.spgrid import spGrid as leafGrid
    padapt = True
except ImportError:
    from dune.grid import yaspGrid as leafGrid
    padapt = False

from dune.grid import cartesianDomain

from dune.fem.space import dglegendrehp as dgspace
from dune.geometry import quadratureRules, quadratureRule
from dune.common import FieldVector
from ufl import SpatialCoordinate, sin, pi
from dune.femdg import createLimiter
import dune.fem as fem

# create one element grid for reference element
domain = cartesianDomain([0], [1], [1])
gridView = leafGrid(domain)

spc = dgspace(gridView, order=6)

x = SpatialCoordinate( spc )
u = 0.5*(sin((x[0]-0.05)*2*pi)*1.1) + 0.5
#u = x[0]*x[0] + x[0] - 0.5

u_h = spc.interpolate( u, name="u_h")
u_h_limited = u_h.copy(name="u_h_limited")

minmax = [1e-12, 1.]

limiter = createLimiter( spc, bounds=[minmax], limiter="scaling" )

def evaluate(df, xval):
    yval = []
    for e in gridView.elements:
        lf = df.localFunction(e)
        for x in xval:
            f = lf.evaluate(x)
            yval.append( f[ 0 ] )
    return np.array(yval)

def evaluateQuad(df):
    yval = []
    xval = []
    for e in gridView.elements:
        lf = df.localFunction(e)
        for p in quadratureRule(e.type, 2*spc.localOrder(e)):
            xval.append(p.position[0])

        xval.sort()

        for x in xval:
            f = lf.evaluate(FieldVector([x]))
            yval.append( f[ 0 ] )
    return np.array(xval), np.array(yval)

def checkLimiter(k, plotting=True):

    markp = lambda E : k

    if padapt:
        # adjust polynomial order
        fem.spaceAdapt( markp, u_h, u_h_limited )

    u_h.interpolate( u )

    x_val = np.linspace(0, 1, 100)

    y_val = evaluate(u_h, x_val)
    y_one = np.ones(len(x_val) )
    y_zero = np.zeros(len(x_val))
    limiter(u_h, u_h_limited)
    y_val = evaluate(u_h_limited, x_val)
    x_quad, y_quad = evaluateQuad(u_h_limited)
    # assure that function values are within bounds
    for y in y_quad:
        assert y >= minmax[0] - 1e-12
        assert y <= minmax[1] + 1e-12

    ma = max(y_val)
    mi = min(y_val)
    #assert ma <= 1.0
    #assert mi >= 0.0
    y_val2 = evaluate(u_h, x_val)
    if plotting:
        plt.plot( x_val, y_one )
        plt.plot( x_val, y_zero)
        plt.plot( x_val, y_val, label=f"u_limited (k={k})" )
        plt.plot( x_quad, y_quad, 'o', label=f"u_limited quad (k={k})" )
        plt.plot( x_val, y_val2, label=f"u_h (k={k})" )
        plt.legend()
        plt.show()

if __name__ == '__main__':
    for k in range(0,spc.order+1):
      checkLimiter(k, plotting=False)
