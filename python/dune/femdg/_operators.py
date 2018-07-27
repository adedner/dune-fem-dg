from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, preprocessorAssert, ConfigurationError

from dune.generator import Constructor, Method
from dune.fem.operator import load

from dune.ufl.tensors import ExprTensor
from dune.ufl.codegen import generateCode
from dune.source.cplusplus import assign, TypeAlias, Declaration, Variable,\
        UnformattedBlock, UnformattedExpression, Struct, return_,\
        SwitchStatement
from dune.source.cplusplus import Method as clsMethod
from dune.source.cplusplus import SourceWriter, ListWriter, StringWriter

from ufl import as_vector, TestFunction,TrialFunction,Coefficient, dx,ds,grad,inner,zero,FacetNormal,dot
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
            dimR = 1
            expr = as_vector([expr])
        t = ExprTensor((dimR, ), [expr[int(i)] for i in range(dimR)])
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
        else:
            predefined = {}
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
        code = [assign(result, construct(cppTye,default) )]
        if returnResult:
            code = [Declaration(result)] + code + [return_(result)]
    return code
def generateMethod(struct,expr, cppType, name,
        returnResult=True,
        defaultReturn='0',
        targs=None, args=None, static=False, const=False, volatile=False,
        predefined={}):
    if not returnResult:
        args = args + [cppType + ' &result']

    if isinstance(expr,dict):
        bndId = Variable('const int', 'bndId')
        code = SwitchStatement(bndId, default=return_(False))
        for id, e in expr.items():
            code.append(id,
                    [generateMethodBody('RangeType', e, False, defaultReturn,
                        predefined), return_(True)])
        code = [code]
    else:
        code = generateMethodBody(cppType, expr, returnResult, defaultReturn,predefined)

    meth = clsMethod(cppType, name,
            code=code,
            args=args,
            targs=targs, static=static, const=const, volatile=volatile)
    struct.append(meth)

#####################################################

# create DG operator + solver
def createFemDGSolver(Model, space):
    import dune.create as create

    u = TrialFunction(space)
    v = TestFunction(space)
    n = FacetNormal(space.cell())
    if hasattr(Model,"F_c"):
        advModel = inner(Model.F_c(u),grad(v))*dx
    else:
        advModel = inner(grad(u-u),grad(v))*dx    # TODO: make a better empty model
    if hasattr(Model,"S_ns"):
        advModel += inner(S_ns(u,grad(u)),v)*dx
    if hasattr(Model,"F_d"):
        diffModel = inner(Model.F_d(u,grad(u)),grad(v))*dx
    else:
        diffModel = inner(grad(u-u),grad(v))*dx   # TODO: make a better empty model
    if hasattr(Model,"S_s"):
        diffModel += inner(S_s(u,grad(u)),v)*dx

    # TODO: needs more general treatment of boundaries
    # advModel  -= inner(dot(Model.F_c(u),n),v)*ds
    # diffModel -= inner(dot(Model.F_c(u),n),v)*ds

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

    includes  = ["dune/fem-dg/solver/dg.hh"]
    includes += space._includes + destinationIncludes
    includes += ["dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"]

    additionalType = 'Additional< typename ' + spaceType + '::FunctionSpaceType >'

    typeName = 'Dune::Fem::DGOperator< ' + destinationType + ', ' + advModelType + ', ' + diffModelType + ', ' + additionalType + ' >'

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
    setTimeStepSize = Method('setTimeStepSize', '&DuneType::setTimeStepSize')

    # add method to obtain time step size
    deltaT = Method('deltaT', '&DuneType::deltaT')

    # extra methods for limiter and time step control
    struct = Struct('Additional', targs=['class FunctionSpace'])
    struct.append(TypeAlias('DomainType','typename FunctionSpace::DomainType'))
    struct.append(TypeAlias('RangeType','typename FunctionSpace::RangeType'))
    struct.append(TypeAlias('JacobianRangeType','typename FunctionSpace::JacobianRangeType'))
    struct.append(TypeAlias('HessianRangeType','typename FunctionSpace::HessianRangeType'))

    arg_n = Variable("const DomainType &", "normal")
    predefined = {n: arg_n}
    maxSpeed = getattr(Model,"maxLambda",None)
    maxSpeed = maxSpeed(u,n)
    generateMethod(struct, maxSpeed,
            'double', 'maxSpeed',
            args=['const Entity &entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    velocity = getattr(Model,"velocity",None)
    velocity = velocity(u)
    generateMethod(struct, velocity,
            'DomainType', 'velocity',
            args=['const Entity &entity', 'const Point &x',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    physical = getattr(Model,"physical",None)
    physical = physical(u)
    generateMethod(struct, physical,
            'double', 'physical',
            args=['const Entity &entity', 'const Point &x',
                  'const RangeType &u'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    w = Coefficient(space)
    predefined.update( {w:Variable("const RangeType &", "w")} )
    jump = getattr(Model,"jump",None)
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
        for id,f in boundaryFluxDict.items(): boundaryFlux.update({id:f(u,n)})
    else:
        boundaryFlux = {}
    generateMethod(struct, boundaryFlux,
            'bool', 'boundaryFlux',
            args=['const int bndId',
                  'const Entity& entity', 'const Point &x',
                  'const DomainType &normal',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)
    boundaryValueDict = getattr(Model,"boundaryValue",None)
    if boundaryValueDict is not None:
        boundaryValue = {}
        for id,f in boundaryValueDict.items(): boundaryValue.update({id:f(u)})
    else:
        boundaryValue = {}
    generateMethod(struct, boundaryValue,
            'bool', 'boundaryValue',
            args=['const int bndId',
                  'const Entity& entity', 'const Point &x',
                  'const RangeType &u',
                  'RangeType &result'],
            targs=['class Entity, class Point'], static=True,
            predefined=predefined)

    writer = SourceWriter(StringWriter())
    writer.emit([struct])

    print("#################################")
    print(writer.writer.getvalue())
    print("#################################")

    return load(includes, typeName, constructor, setTimeStepSize, deltaT,
              preamble=writer.writer.getvalue()).\
                    Operator( space, advModel, diffModel )
