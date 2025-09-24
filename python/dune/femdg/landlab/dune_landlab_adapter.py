import numpy as np
from matplotlib import pyplot as plt
from dune.grid import cartesianDomain
from dune.spgrid import spGrid as gridView
from dune.fem.space import lagrange, finiteVolume, raviartThomas

import time

class DuneLandlabAdapter:
    def __init__(self,grid):

        ## compute lower left and upper right coordinates
        self.lower = [grid.x_of_node[0], grid.y_of_node[0]]
        self.upper = [grid.x_of_node[-1], grid.y_of_node[-1]]

        ## obtain number of cells in each coordinate direction
        self.N = [grid.number_of_node_columns, grid.number_of_node_rows]

        ## compute dx and dy
        self.h = [(self.upper[i] - self.lower[i])/(self.N[i]-1) for i in range(len(self.N))]

        ## add half a cell in each direction
        for i in range(len(self.N)):
            self.lower[ i ] -= 0.5 * self.h[i]
            self.upper[ i ] += 0.5 * self.h[i]

        ## setup DUNE structures
        domain = cartesianDomain(self.lower, self.upper, self.N, overlap=1)
        self.gridView = gridView( domain )

        """
        print(grid.spacing, self.lower, self.upper, self.h, [(u-l)/n for u,l,n in zip(self.upper,self.lower,self.N)])
        linkSpace = raviartThomas( self.gridView, order=0 )
        gradElev = linkSpace.function(name="grad")
        # gradElev.as_numpy[0] = grid.spacing[1]
        gradElev.as_numpy[:] = grid.spacing[1]
        el = self.gridView.elements.__next__()
        vx = np.float64(gradElev.localFunction(el)([0.5,0])[0])
        vy = np.float64(gradElev.localFunction(el)([0.5,0])[1])
        print(gradElev.as_numpy[0],":",vx,vy)
        gradElev[0].plot(block=False)
        gradElev[0].plot(block=False)
        plt.show()
        assert False
        """

        self.grid = grid
        self.linkSpace = dict()    # RT
        self.nodeSpace = dict()    # FV
        self.cellSpace = dict()    # Lagrange

        self._linkIndexMap = None
        self._factor = None

        # setup link identification
        self._setupLinkIndexMap()

    def _setupLinkIndexMap(self):
        #start = time.time()

        factor = [self.grid.spacing[i] for i in range(len(self.N))]

        self._linkIndexMap = np.zeros(self.grid.number_of_links, dtype=int)
        self._factor = np.zeros(self.grid.number_of_links, dtype=float)

        # in yasp first horizontal edges are numbered then vertical edges
        # landlab links are mixed vertical/horizontal
        # In the following l is a counter for landlab data and k is counter for dune data

        # first do vertical edges
        l = 0
        k = self.N[0]*(self.N[1]+1) # skip all horizontal edges
        for i in range(self.N[1]):
            k += 1 # ignore left boundary edge
            for j in range(self.N[0]-1):
                self._linkIndexMap[ l ] = k
                self._factor[l] = factor[0]
                k += 1
                l += 1
            k += 1    # ignore right boundary edge
            l += self.N[0] # skip vertical edges

        # now do horizontal edges
        l = self.N[0]-1 # skip the vertical edges (horizontal links)
        k = self.N[0]      # skip top boundary edges
        for i in range(self.N[1]-1):
            for j in range(self.N[0]):
                self._linkIndexMap[ l ] = k
                self._factor[l] = factor[1]
                k += 1
                l += 1
            l += self.N[0]-1 # skip horizontal edges

        #print(f"Setup took {time.time()-start}s")
        return

    def nodeFct(self,name, dimRange=1):
        if not dimRange in self.nodeSpace:
            self.nodeSpace[dimRange] = finiteVolume( self.gridView, order=0, dimRange=dimRange )
        return self.nodeSpace[ dimRange ].function(name=name)

    def linkFct(self,name, dimRange=1):
        if not dimRange in self.linkSpace:
            self.linkSpace[ dimRange ] = raviartThomas( self.gridView, order=0 )
        return self.linkSpace[ dimRange ].function(name=name)

    def cellFct(self,name, dimRange=1):
        if not dimRange in self.cellSpace:
            self.cellSpace[ dimRange ] = lagrange( self.gridView, order=1, dimRange=dimRange )
        return self.cellSpace[ dimRange ].function(name=name)

    def fromNode(self,d, dh):
        if isinstance(d, (list,tuple)):
            # create function is string was provided
            if isinstance(dh, str):
                dh = self.nodeFct(dh, dimRange=len(d))
            assert len(d) == dh.space.dimRange, f"Expecting a list of size {dh.space.dimRange}"
            stride = dh.space.localBlockSize
            ## d --> dh
            for i in range(dh.space.dimRange):  # update each of the advected scalars
                 dh.as_numpy[i::stride] = d[i][:]
        else:
            if isinstance(dh, str):
                dh = self.nodeFct(dh, dimRange=1)
            assert dh.space.dimRange == 1
            dh.as_numpy[:] = d[:]
        return dh
    ## end fromNode

    def toNode(self, dh, d):
        if isinstance(d, (list,tuple)):
            assert len(d) == dh.space.dimRange, f"Expecting a list of size {dh.space.dimRange}"
            stride = dh.space.localBlockSize
            ## dh --> d
            for i in range(dh.space.dimRange):  # update each of the advected scalars
                 d[i][:] = dh.as_numpy[i::stride]
        else:
            assert dh.space.dimRange == 1
            d[:] = dh.as_numpy[:]
    ## end toNode


    def fromLink(self, links, rtFct):
        if isinstance(links, (list,tuple)):
            assert False, "dimRange > 1 not implemented yet."
        else:
            # create function is string was provided
            if isinstance(rtFct, str):
                rtFct = self.linkFct(rtFct, dimRange=1)
            assert rtFct.space, "This should be a discrete function"
            rtFct.as_numpy[ self._linkIndexMap ] = links[:] * self._factor[:]
        return rtFct

    def toLink(self,rtFct, links):
        assert rtFct.space.dimRange == 1, "dimRange > 1 not implemented yet."
        links[:] = rtFct.as_numpy[ self._linkIndexMap ] / self._factor[:]
