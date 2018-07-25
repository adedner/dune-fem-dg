from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, preprocessorAssert, ConfigurationError

from dune.generator import Constructor, Method
from dune.generator.generator import SimpleGenerator
from dune.generator.generator import SimpleGenerator
from dune.fem.operator import load

generator = SimpleGenerator("Operator", "Dune::FemPy")

# limiter can be ScalingLimiter or FV based limiter with FV type reconstructions for troubled cells
def createLimiter(domainSpace, rangeSpace=None, bounds = [1e-12,1.], limiter='scaling'):
    if rangeSpace is None:
        rangeSpace = domainSpace

    domainSpaceType = domainSpace._typeName
    rangeSpaceType = rangeSpace._typeName

    _, domainFunctionIncludes, domainFunctionType, _, _ = domainSpace.storage
    _, rangeFunctionIncludes, rangeFunctionType, _, _ = rangeSpace.storage

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
