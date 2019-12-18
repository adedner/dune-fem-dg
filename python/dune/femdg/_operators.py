from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, preprocessorAssert, ConfigurationError

from dune.generator import Constructor, Method
from dune.fem.operator import load

from dune.ufl import Constant
from dune.ufl.tensors import ExprTensor
from dune.ufl.codegen import generateCode, generateMethodBody, generateMethod

from dune.source.cplusplus import assign, construct, TypeAlias, Declaration, Variable,\
        UnformattedBlock, UnformattedExpression, Struct, return_,\
        SwitchStatement
from dune.source.cplusplus import Method as clsMethod
from dune.source.cplusplus import SourceWriter, ListWriter, StringWriter

from ufl import SpatialCoordinate,TestFunction,TrialFunction,Coefficient,\
        as_vector, as_matrix,dx,ds,grad,inner,zero,FacetNormal,dot
from ufl.algorithms.analysis import extract_arguments_and_coefficients as coeff
from ufl.differentiation import Grad


# limiter can be ScalingLimiter or FV based limiter with FV type reconstructions for troubled cells
def createLimiter(domainSpace, rangeSpace=None, bounds = [1e-12,1.], limiter='scaling'):
    if rangeSpace is None:
        rangeSpace = domainSpace

    domainSpaceType = domainSpace._typeName
    rangeSpaceType = rangeSpace._typeName

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = domainSpace.storage
    _, rangeFunctionIncludes, rangeFunctionType, _, _, _ = rangeSpace.storage

    includes = ["dune/fem-dg/operator/limiter/limiter.hh"]
    includes += domainSpace._includes + domainFunctionIncludes
    includes += rangeSpace._includes + rangeFunctionIncludes

    typeName = 'Dune::Fem::ScalingLimiter< ' + domainFunctionType + ', ' + rangeFunctionType + ' >'
    # FV type limiter where FV based reconstructions are done
    if limiter is 'fv':
        typeName = 'Dune::Fem::Limiter< ' + domainFunctionType + ', ' + rangeFunctionType + ' >'

    constructor = Constructor(['const '+domainSpaceType + ' &dSpace, const '+rangeSpaceType + ' &rSpace, double lower,double upper'],
                              ['return new ' + typeName + '(dSpace,rSpace,lower,upper );'],
                              ['"dSpace"_a', '"rSpace"_a', '"lower"_a', '"upper"_a',
                               'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])

    # add method activated to inspect limited cells.
    activated = Method('activated', '&'+typeName+'::activated')

    return load(includes, typeName, constructor, activated).Operator( domainSpace, rangeSpace, bounds[0], bounds[1] )

# new method name, only is kept for convenience
def limiter(domainSpace, rangeSpace=None, bounds = [1e-12,1.], limiter='scaling'):
    return createLimiter( domainSpace, rangeSpace, bounds, limiter )

def createOrderRedcution(domainSpace):

    domainSpaceType = domainSpace._typeName

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = domainSpace.storage

    includes = ["dune/fem-dg/operator/common/orderreduction.hh"]
    includes += domainSpace._includes + domainFunctionIncludes

    typeName = 'Dune::Fem::OrderReduction< ' + domainFunctionType + ' >'

    constructor = Constructor(['const '+domainSpaceType + ' &dSpace'],
                              ['return new ' + typeName + '(dSpace);'],
                              ['"dSpace"_a','pybind11::keep_alive< 1, 2 >()'])

    # add method maxRelevantOrder to get max order that is relevant per cell
    # maxRelevantOrder = Method('maxRelevantOrder', '&'+typeName+'::maxRelevantOrder')

    # return load(includes, typeName, constructor, maxRelevantOrder).Operator( domainSpace )
    return load(includes, typeName, constructor).Operator( domainSpace )

#####################################################

#####################################################

# create DG operator + solver (limiter = none,minmod,vanleer,superbee),
# (diffusionScheme = cdg2,br2,ip,nipg,...)
def femDGSolver(Model, space,
        limiter="minmod", diffusionScheme = "cdg2", threading=False,
        initialTime=0.0, parameters={}):
    import dune.create as create

    if limiter is None or limiter is False:
        limiter = "unlimited"

    if diffusionScheme is None or diffusionScheme is False:
        diffusionScheme = "none"

    u = TrialFunction(space)
    v = TestFunction(space)
    n = FacetNormal(space.cell())
    x = SpatialCoordinate(space.cell())
    # this variable needs to be called time in order
    # to work with the implementation in DiffusionModel
    t = Constant(initialTime,"time")
    predefined = {}
    spatial = Variable('const auto', 'y')
    predefined.update( {x: UnformattedExpression('auto', 'entity.geometry().global( Dune::Fem::coordinate( x ) )') })
    arg_n = Variable("const DomainType &", "normal")
    predefined.update( {n: arg_n} )
    arg_t = Variable("const double &", "t")
    predefined.update( {t: arg_t} )

    hasAdvFlux = hasattr(Model,"F_c")
    if hasAdvFlux:
        advModel = inner(Model.F_c(t,x,u),grad(v))*dx
    else:
        advModel = inner(t*grad(u-u),grad(v))*dx    # TODO: make a better empty model
    hasNonStiffSource = hasattr(Model,"S_ns")
    if hasNonStiffSource:
        advModel += inner(as_vector(Model.S_ns(t,x,u,grad(u))),v)*dx

    hasDiffFlux = hasattr(Model,"F_v")
    if hasDiffFlux:
        diffModel = inner(Model.F_v(t,x,u,grad(u)),grad(v))*dx
    else:
        diffModel = inner(t*grad(u-u),grad(v))*dx   # TODO: make a better empty model
    hasStiffSource = hasattr(Model,"S_s")
    if hasStiffSource:
        diffModel += inner(as_vector(Model.S_s(t,x,u,grad(u))),v)*dx

    advModel  = create.model("elliptic",space.grid, advModel)
    diffModel = create.model("elliptic",space.grid, diffModel)

    spaceType = space._typeName

    modelType = "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >"

    advModelType  = modelType
    diffModelType = modelType

    _, destinationIncludes, destinationType, _, _, _ = space.storage

    ###'###############################################
    ### extra methods for limiter and time step control
    struct = Struct('Additional', targs=['class FunctionSpace'])
    struct.append(TypeAlias('DomainType','typename FunctionSpace::DomainType'))
    struct.append(TypeAlias('RangeType','typename FunctionSpace::RangeType'))
    struct.append(TypeAlias('JacobianRangeType','typename FunctionSpace::JacobianRangeType'))
    struct.append(TypeAlias('HessianRangeType','typename FunctionSpace::HessianRangeType'))

    #advFlux2 = getattr(Model,"F_c",None)
    #if advFlux2 is not None:
    #    advFlux2 = advFlux2(t,x,u)
    #generateMethod(struct, advFlux2,
    #        'JacobianRangeType', 'advection',
    #        returnResult=False,
    #        args=['const double &t',
    #              'const Entity &entity',
    #              'const Point &x',
    #              'const RangeType &u'],
    #        targs=['class Entity, class Point'], static=True,
    #        predefined=predefined)

    maxSpeed = getattr(Model,"maxLambda",None)
    if maxSpeed is not None:
        maxSpeed = maxSpeed(t,x,u,n)
    generateMethod(struct, maxSpeed,
            'double', 'maxSpeed',
            args=['const double &t',
                  'const Entity &entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    velocity = getattr(Model,"velocity",None)
    if velocity is not None:
        velocity = velocity(t,x,u)
    generateMethod(struct, velocity,
            'DomainType', 'velocity',
            args=['const double &t',
                  'const Entity &entity', 'const Point &x',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    # TODO: fill in diffusion time step from Model
    diffusionTimeStep = getattr(Model,"maxDiffusion",None)
    if diffusionTimeStep is not None:
        diffusionTimeStep = diffusionTimeStep(t,x,u)
    generateMethod(struct, diffusionTimeStep,
            'double', 'diffusionTimeStep',
            args=['const Entity& entity', 'const Point &x',
                  'const T& circumEstimate', 'const RangeType &u'],
            targs=['class Entity, class Point, class T'], static=True,
            predefined=None)

    #advFlux2 = getattr(Model,"F_c",None)
    #advFlux2 = advFlux2(t,x,u)
    #generateMethod(struct, advFlux2,
    #        'void', 'advection',
    #        args=['const double &t',
    #              'const Entity &entity',
    #              'const Point &x',
    #              'const RangeType &u',
    #              'const JacobianRangeType& du',
    #              'JacobianRangeType& result'],
    #        targs=['class Entity, class Point'], static=True,
    #        predefined=predefined)

    # QUESTION: should `physical` actually depend on x? Perhaps even t?
    physical = getattr(Model,"physical",None)
    if physical is not None:
        physical = physical(u)
    generateMethod(struct, physical,
            'double', 'physical',
            args=['const Entity &entity', 'const Point &x',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    # QUESTION: should `jump` actually depend on x? Perhaps even t?
    w = Coefficient(space)
    predefined.update( {w:Variable("const RangeType &", "w")} )
    jump = getattr(Model,"jump",None)
    if jump is not None:
        jump = jump(u,w)
    generateMethod(struct, jump,
            'double', 'jump',
            args=['const Intersection& it', 'const Point &x',
                  'const RangeType &u',
                  'const RangeType &w'],
            targs=['class Intersection, class Point'], static=True,
            predefined=predefined)

    adjustAverageValue = {}
    generateMethod(struct, adjustAverageValue,
            'void', 'adjustAverageValue',
            args=['const Entity& entity', 'const Point &x',
                  'RangeType &u'],
            targs=['class Entity, class Point'], static=True, evalSwitch=False,
            predefined=None)

    #####################
    ## boundary treatment
    boundaryDict = getattr(Model,"boundary",{})
    boundaryAFlux = {}
    boundaryDFlux = {}
    boundaryValue = {}
    hasBoundaryValue = {}
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
                try:
                    method = f(t,x,u)
                except TypeError:
                    method = f(t,x,u,n,n)
                boundaryValue.update( [ (kk,method) for kk in ids] )
                # add dummy values for hasBoundaryValue method
                bndReturn = True
                hasBoundaryValue.update( [ (kk,bndReturn) for kk in ids] )

    generateMethod(struct, boundaryAFlux,
            'bool', 'boundaryFlux',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)
    generateMethod(struct, boundaryDFlux,
            'bool', 'diffusionBoundaryFlux',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u',
                  'const JacobianRangeType &jac',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)
    generateMethod(struct, boundaryValue,
            'bool', 'boundaryValue',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)
    generateMethod(struct, hasBoundaryValue,
            'bool', 'hasBoundaryValue',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    limiterModifiedDict = getattr(Model,"limitedRange",None)
    if limiterModifiedDict is None:
        limiterModified = {}
        limitedDimRange = "FunctionSpace :: dimRange"
    else:
        limiterModified = {}
        count = 0
        for id,f in limiterModifiedDict.items(): count += 1
        limitedDimRange = str(count)
    generateMethod(struct, limiterModified,
            'void', 'limitedRange',
            args=['LimitedRange& limRange'],
            targs=['class LimitedRange'], static=True, evalSwitch=False,
            predefined=None)

    ##################################
    ## limiter modification size
    struct.append([Declaration(
        Variable("const int", "limitedDimRange = " + limitedDimRange),
        static=True)])
    ##################################
    ## Add 'has*' properties for model
    struct.append([Declaration(
        Variable("const bool", "hasAdvection"), initializer=hasAdvFlux or hasNonStiffSource,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "hasDiffusion"), initializer=hasDiffFlux,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "hasStiffSource"), initializer=hasStiffSource,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "hasNonStiffSource"), initializer=hasNonStiffSource,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "hasFlux"), initializer=hasAdvFlux or hasDiffFlux,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "threading"), initializer=threading,
        static=True)])

    ###################################################
    ## choose details of discretization (i.e. fluxes)
    ## default settings:
    solverId     = "Dune::Fem::Solver::Enum::fem"
    formId       = "Dune::Fem::Formulation::Enum::primal"
    limiterId    = "Dune::Fem::AdvectionLimiter::Enum::limited"
    limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::minmod"
    advFluxId    = "Dune::Fem::AdvectionFlux::Enum::none"
    diffFluxId   = "Dune::Fem::DiffusionFlux::Enum::none"

    if hasDiffFlux:
        diffFluxId = "Dune::Fem::DiffusionFlux::Enum::"+diffusionScheme
    if hasAdvFlux:
        advFluxId  = "Dune::Fem::AdvectionFlux::Enum::llf"
    if limiter.lower() == "unlimited":
        limiterId = "Dune::Fem::AdvectionLimiter::Enum::unlimited"
    # check for different limiter functions (default is minmod)
    if limiter.lower() == "superbee":
        limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::superbee"
    if limiter.lower() == "vanleer":
        limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::vanleer"

    struct.append([Declaration(
        Variable("const Dune::Fem::Solver::Enum", "solverId = " + solverId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::Formulation::Enum", "formId = " + formId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::AdvectionLimiter::Enum", "limiterId = " + limiterId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::AdvectionLimiterFunction::Enum", "limiterFunctionId = " + limiterFctId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::AdvectionFlux::Enum", "advFluxId = " + advFluxId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::DiffusionFlux::Enum", "diffFluxId = " + diffFluxId),
        static=True)])

    writer = SourceWriter(StringWriter())
    writer.emit([struct])

    # print("#################################")
    # print(writer.writer.getvalue())
    # print("#################################")

    ################################################################
    ### Construct DuneType, includes, and extra methods/constructors
    includes  = ["dune/fem-dg/python/operator.hh"]
    includes += ["dune/fem-dg/solver/dg.hh"]
    includes += space._includes + destinationIncludes
    includes += ["dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"]

    additionalType = 'Additional< typename ' + spaceType + '::FunctionSpaceType >'

    typeName = 'Dune::Fem::DGSolver< ' +\
            destinationType + ', ' +\
            advModelType + ', ' + diffModelType + ', ' + additionalType +\
            " >"

    constructor = Constructor(['const '+spaceType + ' &space',
                               'const '+advModelType + ' &advectionModel',
                               'const '+diffModelType + ' &diffusionModel',
                               'const pybind11::dict &parameters'
                              ],
                              ['return new DuneType(space, advectionModel, diffusionModel, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );'],
                              ['"space"_a',
                               '"advectionModel"_a',
                               '"diffusionModel"_a',
                               '"parameters"_a',
                               'pybind11::keep_alive< 1, 2 >()',
                               'pybind11::keep_alive< 1, 3 >()',
                               'pybind11::keep_alive< 1, 4 >()'])

    # add method activated to inspect limited cells.
    applyLimiter = Method('applyLimiter', '''[](
        DuneType &self, typename DuneType::DestinationType &u) {
        self.limit(u); }''' );
    # add method to obtain time step size
    deltaT = Method('deltaT', '&DuneType::deltaT')
    # add method to set a fixed time step
    setTimeStepSize = Method('setTimeStepSize', '&DuneType::setTimeStepSize')
    # add method to solve one step (not requiring u_h_n)
    step = Method('step', '&DuneType::step', extra=['"target"_a'])

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = space.storage
    base = 'Dune::Fem::SpaceOperatorInterface< ' + domainFunctionType + '>'
    return load(includes, typeName,
              constructor, setTimeStepSize, deltaT, applyLimiter, step,
              baseClasses=[base],
              preamble=writer.writer.getvalue()).\
                    Operator( space, advModel, diffModel, parameters=parameters )

# create DG operator + solver (limiter = none,minmod,vanleer,superbee),
# (diffusionScheme = cdg2,br2,ip,nipg,...)
def femDGOperator(Model, space,
        limiter="minmod", diffusionScheme = "cdg2", threading=False,
        initialTime=0.0, parameters={}):
    import dune.create as create

    if limiter is None or limiter is False:
        limiter = "unlimited"

    if diffusionScheme is None or diffusionScheme is False:
        diffusionScheme = "none"

    u = TrialFunction(space)
    v = TestFunction(space)
    n = FacetNormal(space.cell())
    x = SpatialCoordinate(space.cell())
    t = Constant(initialTime,"t")
    predefined = {}
    spatial = Variable('const auto', 'y')
    predefined.update( {x: UnformattedExpression('auto', 'entity.geometry().global( Dune::Fem::coordinate( x ) )') })
    arg_n = Variable("const DomainType &", "normal")
    predefined.update( {n: arg_n} )
    arg_t = Variable("const double &", "t")
    predefined.update( {t: arg_t} )

    hasAdvFlux = hasattr(Model,"F_c")
    if hasAdvFlux:
        advModel = inner(Model.F_c(t,x,u),grad(v))*dx
    else:
        advModel = inner(t*grad(u-u),grad(v))*dx    # TODO: make a better empty model
    hasNonStiffSource = hasattr(Model,"S_ns")
    if hasNonStiffSource:
        advModel += inner(as_vector(Model.S_ns(t,x,u,grad(u))),v)*dx

    hasDiffFlux = hasattr(Model,"F_v")
    if hasDiffFlux:
        diffModel = inner(Model.F_v(t,x,u,grad(u)),grad(v))*dx
    else:
        diffModel = inner(t*grad(u-u),grad(v))*dx   # TODO: make a better empty model
    hasStiffSource = hasattr(Model,"S_s")
    if hasStiffSource:
        diffModel += inner(as_vector(Model.S_s(t,x,u,grad(u))),v)*dx

    advModel  = create.model("elliptic",space.grid, advModel)
    diffModel = create.model("elliptic",space.grid, diffModel)

    spaceType = space._typeName

    modelType = "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >"

    advModelType  = modelType
    diffModelType = modelType

    _, destinationIncludes, destinationType, _, _, _ = space.storage

    ###'###############################################
    ### extra methods for limiter and time step control
    struct = Struct('Additional', targs=['class FunctionSpace'])
    struct.append(TypeAlias('DomainType','typename FunctionSpace::DomainType'))
    struct.append(TypeAlias('RangeType','typename FunctionSpace::RangeType'))
    struct.append(TypeAlias('JacobianRangeType','typename FunctionSpace::JacobianRangeType'))
    struct.append(TypeAlias('HessianRangeType','typename FunctionSpace::HessianRangeType'))

    #advFlux2 = getattr(Model,"F_c",None)
    #if advFlux2 is not None:
    #    advFlux2 = advFlux2(t,x,u)
    #generateMethod(struct, advFlux2,
    #        'JacobianRangeType', 'advection',
    #        returnResult=False,
    #        args=['const double &t',
    #              'const Entity &entity',
    #              'const Point &x',
    #              'const RangeType &u'],
    #        targs=['class Entity, class Point'], static=True,
    #        predefined=predefined)

    maxSpeed = getattr(Model,"maxLambda",None)
    if maxSpeed is not None:
        maxSpeed = maxSpeed(t,x,u,n)
    generateMethod(struct, maxSpeed,
            'double', 'maxSpeed',
            args=['const double &t',
                  'const Entity &entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    velocity = getattr(Model,"velocity",None)
    if velocity is not None:
        velocity = velocity(t,x,u)
    generateMethod(struct, velocity,
            'DomainType', 'velocity',
            args=['const double &t',
                  'const Entity &entity', 'const Point &x',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    # TODO: fill in diffusion time step from Model
    diffusionTimeStep = getattr(Model,"maxDiffusion",None)
    if diffusionTimeStep is not None:
        diffusionTimeStep = diffusionTimeStep(t,x,u)
    generateMethod(struct, diffusionTimeStep,
            'double', 'diffusionTimeStep',
            args=['const Entity& entity', 'const Point &x',
                  'const T& circumEstimate', 'const RangeType &u'],
            targs=['class Entity, class Point, class T'], static=True,
            predefined=None)

    #advFlux2 = getattr(Model,"F_c",None)
    #advFlux2 = advFlux2(t,x,u)
    #generateMethod(struct, advFlux2,
    #        'void', 'advection',
    #        args=['const double &t',
    #              'const Entity &entity',
    #              'const Point &x',
    #              'const RangeType &u',
    #              'const JacobianRangeType& du',
    #              'JacobianRangeType& result'],
    #        targs=['class Entity, class Point'], static=True,
    #        predefined=predefined)

    # QUESTION: should `physical` actually depend on x? Perhaps even t?
    physical = getattr(Model,"physical",None)
    if physical is not None:
        physical = physical(u)
    generateMethod(struct, physical,
            'double', 'physical',
            args=['const Entity &entity', 'const Point &x',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    # QUESTION: should `jump` actually depend on x? Perhaps even t?
    w = Coefficient(space)
    predefined.update( {w:Variable("const RangeType &", "w")} )
    jump = getattr(Model,"jump",None)
    if jump is not None:
        jump = jump(u,w)
    generateMethod(struct, jump,
            'double', 'jump',
            args=['const Intersection& it', 'const Point &x',
                  'const RangeType &u',
                  'const RangeType &w'],
            targs=['class Intersection, class Point'], static=True,
            predefined=predefined)

    adjustAverageValue = {}
    generateMethod(struct, adjustAverageValue,
            'void', 'adjustAverageValue',
            args=['const Entity& entity', 'const Point &x',
                  'RangeType &u'],
            targs=['class Entity, class Point'], static=True, evalSwitch=False,
            predefined=None)

    #####################
    ## boundary treatment
    boundaryDict = getattr(Model,"boundary",{})
    boundaryAFlux = {}
    boundaryDFlux = {}
    boundaryValue = {}
    hasBoundaryValue = {}
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
                try:
                    method = f(t,x,u)
                except TypeError:
                    method = f(t,x,u,n,n)
                boundaryValue.update( [ (kk,method) for kk in ids] )
                # add dummy values for hasBoundaryValue method
                bndReturn = True
                hasBoundaryValue.update( [ (kk,bndReturn) for kk in ids] )

    generateMethod(struct, boundaryAFlux,
            'bool', 'boundaryFlux',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)
    generateMethod(struct, boundaryDFlux,
            'bool', 'diffusionBoundaryFlux',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u',
                  'const JacobianRangeType &jac',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    generateMethod(struct, hasBoundaryValue,
            'bool', 'hasBoundaryValue',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)
    generateMethod(struct, boundaryValue,
            'bool', 'boundaryValue',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    limiterModifiedDict = getattr(Model,"limitedRange",None)
    if limiterModifiedDict is None:
        limiterModified = {}
        limitedDimRange = "FunctionSpace :: dimRange"
    else:
        limiterModified = {}
        count = 0
        for id,f in limiterModifiedDict.items(): count += 1
        limitedDimRange = str(count)
    generateMethod(struct, limiterModified,
            'void', 'limitedRange',
            args=['LimitedRange& limRange'],
            targs=['class LimitedRange'], static=True, evalSwitch=False,
            predefined=None)

    ##################################
    ## limiter modification size
    struct.append([Declaration(
        Variable("const int", "limitedDimRange = " + limitedDimRange),
        static=True)])
    ##################################
    ## Add 'has*' properties for model
    struct.append([Declaration(
        Variable("const bool", "hasAdvection"), initializer=hasAdvFlux or hasNonStiffSource,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "hasDiffusion"), initializer=hasDiffFlux,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "hasStiffSource"), initializer=hasStiffSource,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "hasNonStiffSource"), initializer=hasNonStiffSource,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "hasFlux"), initializer=hasAdvFlux or hasDiffFlux,
        static=True)])
    struct.append([Declaration(
        Variable("const bool", "threading"), initializer=threading,
        static=True)])

    ###################################################
    ## choose details of discretization (i.e. fluxes)
    ## default settings:
    solverId     = "Dune::Fem::Solver::Enum::fem"
    formId       = "Dune::Fem::Formulation::Enum::primal"
    limiterId    = "Dune::Fem::AdvectionLimiter::Enum::limited"
    limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::minmod"
    advFluxId    = "Dune::Fem::AdvectionFlux::Enum::none"
    diffFluxId   = "Dune::Fem::DiffusionFlux::Enum::none"

    if hasDiffFlux:
        diffFluxId = "Dune::Fem::DiffusionFlux::Enum::"+diffusionScheme
    if hasAdvFlux:
        advFluxId  = "Dune::Fem::AdvectionFlux::Enum::llf"
    if limiter.lower() == "unlimited":
        limiterId = "Dune::Fem::AdvectionLimiter::Enum::unlimited"
    # check for different limiter functions (default is minmod)
    if limiter.lower() == "superbee":
        limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::superbee"
    if limiter.lower() == "vanleer":
        limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::vanleer"

    struct.append([Declaration(
        Variable("const Dune::Fem::Solver::Enum", "solverId = " + solverId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::Formulation::Enum", "formId = " + formId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::AdvectionLimiter::Enum", "limiterId = " + limiterId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::AdvectionLimiterFunction::Enum", "limiterFunctionId = " + limiterFctId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::AdvectionFlux::Enum", "advFluxId = " + advFluxId),
        static=True)])
    struct.append([Declaration(
        Variable("const Dune::Fem::DiffusionFlux::Enum", "diffFluxId = " + diffFluxId),
        static=True)])

    writer = SourceWriter(StringWriter())
    writer.emit([struct])

    # print("#################################")
    # print(writer.writer.getvalue())
    # print("#################################")

    ################################################################
    ### Construct DuneType, includes, and extra methods/constructors
    includes  = ["dune/fem-dg/python/operator.hh"]
    includes += ["dune/fem-dg/operator/dg/dgpyoperator.hh"]
    includes += space._includes + destinationIncludes
    includes += ["dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"]

    additionalType = 'Additional< typename ' + spaceType + '::FunctionSpaceType >'

    typeName = 'Dune::Fem::DGOperator< ' +\
            destinationType + ', ' +\
            advModelType + ', ' + diffModelType + ', ' + additionalType +\
            " >"

    constructor = Constructor(['const '+spaceType + ' &space',
                               'const '+advModelType + ' &advectionModel',
                               'const '+diffModelType + ' &diffusionModel',
                               'const pybind11::dict &parameters'
                              ],
                              ['return new DuneType(space, advectionModel, diffusionModel, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );'],
                              ['"space"_a',
                               '"advectionModel"_a',
                               '"diffusionModel"_a',
                               '"parameters"_a',
                               'pybind11::keep_alive< 1, 2 >()',
                               'pybind11::keep_alive< 1, 3 >()',
                               'pybind11::keep_alive< 1, 4 >()'])

    # add method activated to inspect limited cells.
    applyLimiter = Method('applyLimiter', '''[](
        DuneType &self, typename DuneType::DestinationType &u) {
        self.limit(u); }''' );

    # add method to extract explicit and implicit part of the operator
    #explOp  = Method('explicitOperator', '&DuneType::explicitOperator')
    #implOp  = Method('implicitOperator', '&DuneType::implicitOperator')

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = space.storage
    base = 'Dune::Fem::SpaceOperatorInterface< ' + domainFunctionType + '>'
    return load(includes, typeName,
                baseClasses = [base],
                preamble=writer.writer.getvalue()).\
                Operator( space, advModel, diffModel, parameters=parameters )

# RungeKutta solvers
def rungeKuttaSolver( fullOperator, imex='EX', butchertable=None, parameters={} ):

    includes = ["dune/fem-dg/solver/rungekuttasolver.hh", "dune/fem-dg/misc/algorithmcreatorselector.hh"]
    includes += fullOperator._includes

    space = fullOperator.domainSpace
    spaceType = space._typeName

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = space.storage

    baseOperatorType = 'Dune::Fem::SpaceOperatorInterface< ' + domainFunctionType + '>'
    fullOperatorType = baseOperatorType
    explOperatorType = baseOperatorType
    implOperatorType = baseOperatorType

    typeName = 'Dune::Fem::SimpleRungeKuttaSolver< ' + domainFunctionType + '>'

    imexId = 0 # == 'EX'
    if imex == 'IM':
        imexId = 1
    elif imex == 'IMEX':
        imexId = 2

    # TODO: move this to header file in dune/fem-dg/python
    # constructor = Constructor([operatorType + ' &op'],
    #                           ['return new ' + typeName + '(op);'],
    #                           ['"op"_a','pybind11::keep_alive< 1, 2 >()' ])
    constructor = Constructor([fullOperatorType + ' &op',
                               explOperatorType + ' &explOp',
                               implOperatorType + ' &implOp',
                               'const int imexId',
                               'const pybind11::dict &parameters'],
                              ['return new ' + typeName + '(op, explOp, implOp, imexId, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ));'],
                              ['"op"_a', '"explOp"_a', '"implOp"_a', '"imexId"_a', '"parameters"_a',
                               'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()','pybind11::keep_alive< 1, 4 >()' ])

    solve = Method('step', '''[]( DuneType &self, typename DuneType::DestinationType &u) { self.solve(u); }''' );
    # add method to solve one step (not requiring u_h_n)
    step = Method('step', '&DuneType::step', extra=['"target"_a'])
    setTimeStepSize = Method('setTimeStepSize', '&DuneType::setTimeStepSize')
    deltaT = Method('deltaT', '&DuneType::deltaT')

    return load(includes, typeName, constructor, step, solve, setTimeStepSize, deltaT).Operator(
            fullOperator.fullOperator,
            fullOperator.explicitOperator,
            fullOperator.implicitOperator,
            imexId,
            parameters=parameters)
