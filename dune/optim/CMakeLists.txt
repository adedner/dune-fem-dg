add_subdirectory( common )
add_subdirectory( constraint )
add_subdirectory( qp )
add_subdirectory( solver )
add_subdirectory( std )

set( HEADERS
  activeindexmapper.hh
  lp.hh
  minimize.hh
  qp.hh
)

install( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/optim )
