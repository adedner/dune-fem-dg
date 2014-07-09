set(miscdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/misc)
set(misc_HEADERS  diagnostics.hh cons2prim.hh crs.hh streams.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${misc_HEADERS} DESTINATION ${miscdir})