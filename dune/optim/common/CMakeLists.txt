set( HEADERS
  axisalignedreferencefaces.hh
  densesubvector.hh
  fieldmatrix.hh
  getreference.hh
  integerrange.hh
  matrix.hh
  smallobject.hh
)

install( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/optim/common )
