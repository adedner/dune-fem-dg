set(solverdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/solver)
set(solver_HEADERS  smartodesolver.hh  rungekuttasolver.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${solver_HEADERS} DESTINATION ${solverdir})