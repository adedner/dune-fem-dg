from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, preprocessorAssert, ConfigurationError

from dune.generator import Constructor, Method
from dune.fem.operator import load

from dune.ufl import NamedConstant
from dune.ufl.tensors import ExprTensor
from dune.ufl.codegen import generateCode
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

def generateMethodBody(cppType, expr, returnResult, default, predefined):
    if expr is not None:
        try:
            dimR = expr.ufl_shape[0]
        except:
            if isinstance(expr,list) or isinstance(expr,tuple):
                expr = as_vector(expr)
            else:
                expr = as_vector([expr])
            dimR = expr.ufl_shape[0]
        # t = ExprTensor((dimR, ), [expr[int(i)] for i in range(dimR)])
        t = ExprTensor(expr.ufl_shape, [expr[int(i)] for i in range(dimR)])
        expression = [expr[i] for i in t.keys()]
        u = coeff(expr)[0]
        if u != []:
            u = u[0]
            du = Grad(u)
            d2u = Grad(du)
            arg_u = Variable("const RangeType &", "u")
            arg_du = Variable("const JacobianRangeType &", "du")
            arg_d2u = Variable("const HessianRangeType &", "d2u")
            predefined.update( {u: arg_u, du: arg_du, d2u: arg_d2u} )
        code, results = generateCode(predefined, expression, tempVars=False)
        result = Variable(cppType, 'result')
        if cppType == 'double':
            code = code + [assign(result, results[0])]
        else:
            code = code + [assign(result[i], r) for i, r in zip(t.keys(), results)]
        if returnResult:
            code = [Declaration(result)] + code + [return_(result)]
    else:
        result = Variable(cppType, 'result')
        code = [assign(result, construct(cppType,default) )]
        if returnResult:
            code = [Declaration(result)] + code + [return_(result)]
    return code
def generateMethod(struct,expr, cppType, name,
        returnResult=True,
        defaultReturn='0',
        targs=None, args=None, static=False, const=False, volatile=False,
        evalSwitch=True,
        predefined={}):
    if not returnResult:
        args = args + [cppType + ' &result']
        returnType = 'void'
    else:
        returnType = cppType

    if isinstance(expr,dict):
        if evalSwitch:
            bndId = Variable('const int', 'bndId')
            code = SwitchStatement(bndId, default=return_(False))
            for id, e in expr.items():
                code.append(id,
                        [generateMethodBody('RangeType', e, False, defaultReturn,
                            predefined), return_(True)])
        else:
            code = UnformattedBlock()
        code = [code]
    else:
        code = generateMethodBody(cppType, expr, returnResult, defaultReturn, predefined)

    meth = clsMethod(returnType, name,
            code=code,
            args=args,
            targs=targs, static=static, const=const, volatile=volatile)
    struct.append(meth)

#####################################################

# create DG operator + solver
def createFemDGSolver(Model, space,
        limiter="default", diffusionScheme = "cdg2"):
    import dune.create as create

    u = TrialFunction(space)
    v = TestFunction(space)
    n = FacetNormal(space.cell())
    x = SpatialCoordinate(space.cell())
    t = NamedConstant(space,"t")
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
        advModel += inner(as_vector(S_ns(t,x,u,grad(u))),v)*dx

    hasDiffFlux = hasattr(Model,"F_v")
    if hasDiffFlux:
        diffModel = inner(Model.F_v(t,x,u,grad(u)),grad(v))*dx
    else:
        diffModel = inner(t*grad(u-u),grad(v))*dx   # TODO: make a better empty model
    hasStiffSource = hasattr(Model,"S_s")
    if hasStiffSource:
        diffModel += inner(as_vector(S_s(t,x,u,grad(u))),v)*dx

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

    advFlux2 = getattr(Model,"F_c",None)
    if advFlux2 is not None:
        advFlux2 = advFlux2(t,x,u)
    generateMethod(struct, advFlux2,
            'JacobianRangeType', 'advection',
            returnResult=False,
            args=['const double &t',
                  'const Entity &entity',
                  'const Point &x',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)


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

    boundaryFluxDict = getattr(Model,"boundaryFlux",None)
    if boundaryFluxDict is not None:
        boundaryFlux = {}
        for id,f in boundaryFluxDict.items(): boundaryFlux.update({id:f(t,x,u,n)})
    else:
        boundaryFlux = {}
    generateMethod(struct, boundaryFlux,
            'bool', 'boundaryFlux',
            args=['const int bndId',
                  'const double &t',
                  'const Entity& entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)
    diffusionBoundaryFluxDict = getattr(Model,"diffusionBoundaryFlux",None)
    if diffusionBoundaryFluxDict is not None:
        diffusionBoundaryFlux = {}
        for id,f in diffusionBoundaryFluxDict.items(): diffusionBoundaryFlux.update({id:f(t,x,u,n)})
    else:
        diffusionBoundaryFlux = {}
    generateMethod(struct, diffusionBoundaryFlux,
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
    boundaryValueDict = getattr(Model,"boundaryValue",None)
    if boundaryValueDict is not None:
        boundaryValue = {}
        for id,f in boundaryValueDict.items(): boundaryValue.update({id:f(t,x,u)})
    else:
        boundaryValue = {}
    generateMethod(struct, boundaryValue,
            'bool', 'boundaryValue',
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

    ###################################################
    ## choose details of discretization (i.e. fluxes)
    ## default settings:
    solverId   = "Dune::Fem::Solver::Enum::fem"
    formId     = "Dune::Fem::Formulation::Enum::primal"
    limiterId  = "Dune::Fem::AdvectionLimiter::Enum::limited"
    advFluxId  = "Dune::Fem::AdvectionFlux::Enum::none"
    diffFluxId = "Dune::Fem::DiffusionFlux::Enum::none"

    if hasDiffFlux:
        diffFluxId = "Dune::Fem::DiffusionFlux::Enum::"+diffusionScheme
    if hasAdvFlux:
        advFluxId  = "Dune::Fem::AdvectionFlux::Enum::llf"
    if limiter == None or limiter == False or limiter.lower() == "unlimiter":
        limiterId = "Dune::Fem::AdvectionLimiter::Enum::unlimited"

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
    includes  = ["dune/fem-dg/solver/dg.hh"]
    includes += space._includes + destinationIncludes
    includes += ["dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"]

    additionalType = 'Additional< typename ' + spaceType + '::FunctionSpaceType >'

    typeName = 'Dune::Fem::DGOperator< ' +\
            destinationType + ', ' +\
            advModelType + ', ' + diffModelType + ', ' + additionalType +\
            " >"

    constructor = Constructor(['const '+spaceType + ' &space',
                               'const '+advModelType + ' &advectionModel',
                               'const '+diffModelType + ' &diffusionModel'
                              ],
                              ['return new DuneType(space, advectionModel, diffusionModel);'],
                              ['"space"_a',
                               '"advectionModel"_a',
                               '"diffusionModel"_a',
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
    # add method to solve (not requiring u_h_n)
    solve = Method('solve', '&DuneType::solve', extra=['"target"_a'])

    return load(includes, typeName, constructor, setTimeStepSize, deltaT, applyLimiter, solve,
              preamble=writer.writer.getvalue()).\
                    Operator( space, advModel, diffModel )
