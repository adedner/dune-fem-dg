set(operatordir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/operator/limiter)
set(operator_HEADERS  limiter.hh limitpass.hh limiterdiscretemodel.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${operator_HEADERS} DESTINATION ${operatordir})