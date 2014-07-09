set(passdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/pass)
set(pass_HEADERS  assembleddiffusionpass.hh  dgmodelcaller.hh  dgmasspass.hh
               dgpass.hh  ellipticmodelcaller.hh  threadhandle.hh  threadpass.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${pass_HEADERS} DESTINATION ${passdir})