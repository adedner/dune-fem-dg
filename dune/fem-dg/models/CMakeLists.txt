set(modelsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/models)
set(models_HEADERS  defaultmodel.hh defaultprobleminterfaces.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${models_HEADERS} DESTINATION ${modelsdir})