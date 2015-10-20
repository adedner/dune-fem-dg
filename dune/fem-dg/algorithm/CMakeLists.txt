set(algorithmsrcdir  ${top_srcdir}/dune/fem/algorithm)
set(algorithmsrc_HEADERS  algorithm.hh  base.hh  evolution.hh  subevolution.hh  monitor.hh substeadystate.hh steadystate.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${algorithmsrc_HEADERS} DESTINATION ${algorithmsrcdir})
