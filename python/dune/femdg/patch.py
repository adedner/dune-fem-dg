from __future__ import division, print_function, unicode_literals

from dune.ufl.codegen import generateMethod
from ufl import grad, TrialFunction, SpatialCoordinate, FacetNormal, Coefficient, replace, diff, as_vector
from ufl.core.expr import Expr
from dune.source.cplusplus import Variable, UnformattedExpression, AccessModifier
from ufl.algorithms import expand_compounds, expand_derivatives, expand_indices, expand_derivatives

def codeFemDg(self):
    code = self._code()
    code.append(AccessModifier("public"))

    space = self._space
    Model = self._Model
    t = self._t

    u = TrialFunction(space)
    # v = TestFunction(space)
    n = FacetNormal(space.cell())
    x = SpatialCoordinate(space.cell())

    predefined = {}
    spatial = Variable('const auto', 'y')
    predefined.update( {x: UnformattedExpression('auto', 'entity.geometry().global( Dune::Fem::coordinate( x ) )') })
    arg_n = Variable("const DDomainType &", "normal")
    predefined.update( {n: arg_n} )
    arg_t = Variable("const double &", "t")
    predefined.update( {t: arg_t} )

    maxSpeed = getattr(Model,"maxLambda",None)
    if maxSpeed is not None:
        maxSpeed = maxSpeed(t,x,u,n)
    generateMethod(code, maxSpeed,
            'double', 'maxSpeed',
            args=['const double &t',
                  'const Entity &entity', 'const Point &x',
                  'const DDomainType &normal',
                  'const DRangeType &u'],
            targs=['class Entity, class Point'], const=True, static=False,
            predefined=predefined)

    velocity = getattr(Model,"velocity",None)
    if velocity is not None:
        velocity = velocity(t,x,u)
    generateMethod(code, velocity,
            'DDomainType', 'velocity',
            args=['const double &t',
                  'const Entity &entity', 'const Point &x',
                  'const DRangeType &u'],
            targs=['class Entity, class Point'], const=True, static=False,
            predefined=predefined)

    # TODO: fill in diffusion time step from Model
    diffusionTimeStep = getattr(Model,"maxDiffusion",None)
    if diffusionTimeStep is not None:
        diffusionTimeStep = diffusionTimeStep(t,x,u)
    generateMethod(code, diffusionTimeStep,
            'double', 'diffusionTimeStep',
            args=['const Entity& entity', 'const Point &x',
                  'const T& circumEstimate', 'const DRangeType &u'],
            targs=['class Entity, class Point, class T'], const=True, static=False,
            predefined=None)

    # QUESTION: should `physical` actually depend on x? Perhaps even t?
    physical = getattr(Model,"physical",None)
    if physical is not None:
        physical = physical(u)
    generateMethod(code, physical,
            'double', 'physical',
            args=['const Entity &entity', 'const Point &x',
                  'const DRangeType &u'],
            targs=['class Entity, class Point'], const=True, static=False,
            predefined=predefined)

    # QUESTION: should `jump` actually depend on x? Perhaps even t?
    w = Coefficient(space)
    predefined.update( {w:Variable("const DRangeType &", "w")} )
    jump = getattr(Model,"jump",None)
    if jump is not None:
        jump = jump(u,w)
    generateMethod(code, jump,
            'double', 'jump',
            args=['const Intersection& it', 'const Point &x',
                  'const DRangeType &u',
                  'const DRangeType &w'],
            targs=['class Intersection, class Point'], const=True, static=False,
            predefined=predefined)

    adjustAverageValue = {}
    generateMethod(code, adjustAverageValue,
            'void', 'adjustAverageValue',
            args=['const Entity& entity', 'const Point &x',
                  'DRangeType &u'],
            targs=['class Entity, class Point'], const=True, static=False, evalSwitch=False,
            predefined=None)

    #####################
    ## boundary treatment
    hasAdvFlux = hasattr(Model,"F_c")
    hasDiffFlux = hasattr(Model,"F_v")
    boundaryDict = getattr(Model,"boundary",{})
    boundaryAFlux = {}
    boundaryDFlux = {}
    boundaryValue = {}
    for k,f in boundaryDict.items():
        # collect all ids (could be list or range)
        ids = []
        try:
            for kk in k:
                ids += [kk]
        except TypeError:
            ids += [k]
        # TODO: check that id is not used more then once
        # figure out what type of boundary condition is used
        if isinstance(f,tuple) or isinstance(f,list):
            assert hasAdvFlux and hasDiffFlux, "two boundary fluxes provided for id "+str(k)+" but only one bulk flux given"
            method = [f[0](t,x,u,n), f[1](t,x,u,grad(u),n)]
            boundaryAFlux.update( [ (kk,method[0]) for kk in ids] )
            boundaryDFlux.update( [ (kk,method[1]) for kk in ids] )
        else:
            try:
                method = f(t,x,u,n)
                if hasAdvFlux and not hasDiffFlux:
                    boundaryAFlux.update( [ (kk,method) for kk in ids] )
                elif not hasAdvFlux and hasDiffFlux:
                    boundaryDFlux.update( [ (kk,method) for kk in ids] )
                else:
                    assert not (hasAdvFlux and hasDiffFlux), "one boundary fluxes provided for id "+str(k)+" but two bulk flux given"
            except TypeError:
                method = f(t,x,u)
                boundaryValue.update( [ (kk,method) for kk in ids] )

    generateMethod(code, boundaryAFlux,
            'bool', 'boundaryFlux',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DDomainType &normal',
                  'const DRangeType &u',
                  'RRangeType &result'],
            targs=['class Entity, class Point'], const=True, static=False,
            predefined=predefined)
    generateMethod(code, boundaryDFlux,
            'bool', 'diffusionBoundaryFlux',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DDomainType &normal',
                  'const DRangeType &u',
                  'const DJacobianRangeType &jac',
                  'RRangeType &result'],
            targs=['class Entity, class Point'], const=True, static=False,
            predefined=predefined)
    generateMethod(code, boundaryValue,
            'bool', 'boundaryValue',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DRangeType &u',
                  'RRangeType &result'],
            targs=['class Entity, class Point'], const=True, static=False,
            predefined=predefined)

    limiterModifiedDict = getattr(Model,"limitedRange",None)
    if limiterModifiedDict is None:
        limiterModified = {}
        limitedDimRange = "DFunctionSpace :: dimRange"
    else:
        limiterModified = {}
        count = 0
        for id,f in limiterModifiedDict.items(): count += 1
        limitedDimRange = str(count)
    generateMethod(code, limiterModified,
            'void', 'limitedRange',
            args=['LimitedRange& limRange'],
            targs=['class LimitedRange'], const=True, static=False, evalSwitch=False,
            predefined=None)

    return code

def transform(Model,space,t):
    def transform_(model):
        if model.baseName == "modelFemDg":
            return
        model._code = model.code
        model.code  = lambda *args,**kwargs: codeFemDg(*args,**kwargs)
        model.baseName = "modelFemDg"
        model._Model = Model
        model._space = space
        model._t = t
        # model.modelWrapper = "DiffusionModelWrapper< Model >"
    return [transform_,None] # add ufl forms in second argument to extract extra coeffs
