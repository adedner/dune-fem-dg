set(SUBDIRS  assemble main misc models operator pass solver stepper test)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
set(femdgdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg)
set(femdg_HEADERS  femdgtest.hh)
install(FILES ${femdg_HEADERS} DESTINATION ${femdgdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
