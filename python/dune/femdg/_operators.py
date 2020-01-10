from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, preprocessorAssert, ConfigurationError

from dune.generator import Constructor, Method
from dune.fem.operator import load
from dune.fem import parameter as parameterReader

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
## fem-dg Operator
#####################################################
from dune.femdg.patch import transform
# create DG operator + solver (limiter = none,minmod,vanleer,superbee),
# (diffusionScheme = cdg2,br2,ip,nipg,...)
def femDGOperator(Model, space,
        limiter="minmod", diffusionScheme = "cdg2", threading=False,
        initialTime=0.0, parameters={}):
    virtualize = False
    import dune.create as create

    if limiter is None or limiter is False:
        limiter = "unlimited"

    if limiter.lower() == "default":
        limiter = "minmod"

    # TODO: does this make sense - if there is no diffusion then it doesn't
    # matter and with diffusion using 'none' seems a bad idea?
    if diffusionScheme is None or diffusionScheme is False:
        diffusionScheme = "none"

    u = TrialFunction(space)
    v = TestFunction(space)
    n = FacetNormal(space.cell())
    x = SpatialCoordinate(space.cell())
    t = Constant(initialTime,"time")

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

    advModel  = create.model("elliptic",space.grid, advModel,
                      modelPatch=transform(Model,space,t,"Adv"),
                      virtualize=virtualize)
    diffModel = create.model("elliptic",space.grid, diffModel,
                      modelPatch=transform(Model,space,t,"Diff"),
                      virtualize=virtualize)

    spaceType = space._typeName

    if virtualize:
        modelType = "DiffusionModel< " +\
              "typename " + spaceType + "::GridPartType, " +\
              spaceType + "::dimRange, " +\
              spaceType + "::dimRange, " +\
              "typename " + spaceType + "::RangeFieldType >"
        advModelType  = modelType
        diffModelType = modelType
    else:
        advModelType  = advModel._typeName # modelType
        diffModelType = diffModel._typeName # modelType

    _, destinationIncludes, destinationType, _, _, _ = space.storage

    ###'###############################################
    ### extra methods for limiter and time step control
    struct = Struct('Additional', targs=['class FunctionSpace'])
    struct.append(TypeAlias('DomainType','typename FunctionSpace::DomainType'))
    struct.append(TypeAlias('RangeType','typename FunctionSpace::RangeType'))
    struct.append(TypeAlias('JacobianRangeType','typename FunctionSpace::JacobianRangeType'))
    struct.append(TypeAlias('HessianRangeType','typename FunctionSpace::HessianRangeType'))

    ##################################
    ## limiter modification size
    limiterModifiedDict = getattr(Model,"limitedRange",None)
    if limiterModifiedDict is None:
        limiterModified = {}
        limitedDimRange = "FunctionSpace :: dimRange"
    else:
        limiterModified = {}
        count = 0
        for id,f in limiterModifiedDict.items(): count += 1
        limitedDimRange = str(count)
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
        # if dgadvectionflux.method has been selected, then use general flux,
        # otherwise default to LLF flux
        key = 'dgadvectionflux.method'
        if key in parameters.keys():
            value = parameters["dgadvectionflux.method"]
            # set parameter in dune-fem parameter container
            parameterReader.append( { key: value } )
            if value.upper().find( 'LLF' ) >= 0:
                advFluxId  = "Dune::Fem::AdvectionFlux::Enum::llf"
            else:
                if value.upper().find( 'EULER' ) >= 0:
                    advFluxId  = "Dune::Fem::AdvectionFlux::Enum::euler_general"
                else:
                    advFluxId  = "Dune::Fem::AdvectionFlux::Enum::general"
        else:
            advFluxId  = "Dune::Fem::AdvectionFlux::Enum::llf"

    if limiter.lower() == "unlimited":
        limiterId = "Dune::Fem::AdvectionLimiter::Enum::unlimited"
    elif limiter.lower() == "scaling":
        limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::none"
        limiterId = "Dune::Fem::AdvectionLimiter::Enum::scalinglimited"
    # check for different limiter functions (default is minmod)
    elif limiter.lower() == "superbee":
        limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::superbee"
    elif limiter.lower() == "vanleer":
        limiterFctId = "Dune::Fem::AdvectionLimiterFunction::Enum::vanleer"
    elif limiter.lower() != "minmod":
        raise ValueError("limiter "+limiter+" not recognised")

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
    includes += advModel._includes + diffModel._includes

    additionalType = 'Additional< typename ' + spaceType + '::FunctionSpaceType >'

    typeName = 'Dune::Fem::DGOperator< ' +\
            destinationType + ', ' +\
            advModelType + ', ' + diffModelType + ', ' + additionalType +\
            " >"

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = space.storage
    base = 'Dune::Fem::SpaceOperatorInterface< ' + domainFunctionType + '>'
    if parameters is not None:
        op = load(includes, typeName,
                 baseClasses = [base],
                 preamble=writer.writer.getvalue()).\
                 Operator( space, advModel, diffModel, parameters=parameters )
    else:
        op = load(includes, typeName,
                 baseClasses = [base],
                 preamble=writer.writer.getvalue()).\
                 Operator( space, advModel, diffModel )
    op._t = t
    op.time = t.value
    op.models = [advModel,diffModel]
    op.space = space
    def setTime(self,time):
        self._t.value = time
        self.time = time
        self._setTime(self.time)
    op.setTime = setTime.__get__(op)
    # def addToTime(self,dt):
    #     self.setTime(self,self.time+dt)
    # op.addToTime = addToTime.__get__(op)
    def stepTime(self,c,dt):
        self.setTime(self.time+c*dt)
    op.stepTime  = stepTime.__get__(op)
    op._hasAdvFlux = hasAdvFlux
    op._hasDiffFlux = hasDiffFlux
    return op

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
    constructor = Constructor([fullOperatorType + ' &op',
                               explOperatorType + ' &explOp',
                               implOperatorType + ' &implOp',
                               'const int imexId',
                               'const pybind11::dict &parameters'],
                              ['return new ' + typeName + '(op, explOp, implOp, imexId, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ));'],
                              ['"op"_a', '"explOp"_a', '"implOp"_a', '"imexId"_a', '"parameters"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()','pybind11::keep_alive< 1, 4 >()','pybind11::keep_alive< 1, 6>()' ])

    solve = Method('solve', '''[]( DuneType &self, typename DuneType::DestinationType &u) { self.solve(u); }''' );
    setTimeStepSize = Method('setTimeStepSize', '&DuneType::setTimeStepSize')
    deltaT = Method('deltaT', '&DuneType::deltaT')

    return load(includes, typeName, constructor, solve, setTimeStepSize, deltaT).Operator(
            fullOperator.fullOperator,
            fullOperator.explicitOperator,
            fullOperator.implicitOperator,
            imexId,
            parameters=parameters)
