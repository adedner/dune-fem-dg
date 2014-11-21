set(SUBDIRS  dataio advdiff euler navierstokes poisson)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
set(testdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/test)
set(test_HEADERS  )
install(FILES ${test_HEADERS} DESTINATION ${testdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
