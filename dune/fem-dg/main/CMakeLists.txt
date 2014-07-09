set(maindir  ${top_srcdir}/dune/fem-dg/main)
set(main_HEADERS  main.hh codegen.hh  codegen2.hh caching.hh caching2.hh
  vectorialbasefunctionset.hh default.hh
  main_0.cc main_1.cc  main_2.cc  main_3.cc  main_4.cc  main.cc main_pol.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${main_HEADERS} DESTINATION ${maindir})