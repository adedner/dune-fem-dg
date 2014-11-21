set(assembledir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/assemble)
set(assemble_HEADERS  primalmatrix.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${assemble_HEADERS} DESTINATION ${assembledir})