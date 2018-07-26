import time
from dune.grid import structuredGrid
from dune.fem import parameter
import dune.create as create
from dune.models.elliptic.formfiles import loadModels
from llf import NumFlux
from ufl import *
from dune.femdg import createFemDGSolver

exec(open("euler.ufl").read())

parameter.append({"fem.verboserank": -1})

grid = structuredGrid([-1, 0], [1, 0.1], [200, 10])
dimR      = 4
UL = as_vector( Model.toCons([1,0,0,1]) )
UR = as_vector( Model.toCons([0.125,0,0,0.1]) )
x = SpatialCoordinate(space.cell())
initial = conditional(x[0]<0,UL,UR)
t = 0
saveStep = 0.01
saveTime = saveStep
count = 0
endTime = 0.4

def useGalerkinOp():
    # a full fv implementation using UFL and the galerkin operator
    # some sign error or something else leading to wrong results
    # still needs fixing
    spaceName = "finitevolume"
    space = create.space(spaceName, grid, dimrange=dimR)

    n = FacetNormal(space.cell())

    u_h   = space.interpolate(initial, name='u_h')
    u_h_n = u_h.copy(name="previous")

    fullModel = inner( Model.F_c(u), grad(v) ) * dx -\
                inner( NumFlux(Model,u,n,dt), v('+')-v('-')) * dS -\
                inner( Model.F_c(u)*n, v) * ds
    operator = create.operator("galerkin",
                 inner(u-u_h_n,v)*dx == dt*replace(fullModel,{u:u_h_n}),
                 space )

    operator.model.dt = 1e-5

    start = time.time()
    grid.writeVTK('sod', pointdata=[u_h], number=count)
    while t < 0.4:
        u_h_n.assign(u_h)
        operator.solve(target=u_h)
        t += operator.model.dt
        if t > saveTime:
            count += 1
            grid.writeVTK('sod', pointdata=[u_h], number=count)
            saveTime += saveStep
    grid.writeVTK('sod', pointdata=[u_h], number=count)
    print("time loop:",time.time()-start)

def useODESolver():
    spaceName = "dgonb"
    polOrder = 0
    space = create.space(spaceName, grid, order=polOrder, dimrange=dimR)
    u_h   = space.interpolate(initial, name='u_h')
    u_h_n = u_h.copy(name="previous")
    operator = createFemDGSolver( Compressible2DEuler, space )

    start = time.time()
    grid.writeVTK('sod', pointdata=[u_h], number=count)
    while t < 0.4:
        u_h_n.assign(u_h)
        operator(u_h_n, u_h)
        t += dt
        if t > saveTime:
            count += 1
            grid.writeVTK('sod', pointdata=[u_h], number=count)
            saveTime += saveStep
    grid.writeVTK('sod', pointdata=[u_h], number=count)
    print("time loop:",time.time()-start)

useODESolver()
