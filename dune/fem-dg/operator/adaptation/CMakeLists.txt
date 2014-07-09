set(adaptationdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/operator/adaptation)
set(adaptation_HEADERS  adaptation.cc  adaptation.hh
                                                                                 estimatorbase.hh estimator.hh utility.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${adaptation_HEADERS} DESTINATION ${adaptationdir})