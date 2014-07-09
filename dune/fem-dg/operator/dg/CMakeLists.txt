set(dgdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/operator/dg)
set(dg_HEADERS  discretemodelcommon.hh  fluxoperator.hh  primaldiscretemodel.hh
dgoperatorchoice.hh  fluxdiscretemodel.hh    operatorbase.hh  primaloperator.hh
passtraits.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${dg_HEADERS} DESTINATION ${dgdir})