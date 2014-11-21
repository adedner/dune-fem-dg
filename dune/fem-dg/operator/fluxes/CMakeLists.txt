set(fluxesdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/operator/fluxes)
set(fluxes_HEADERS  averageflux.hh  diffusionflux.hh  ldgflux.hh
mhd_fluxes.hh  rotator.hh dgprimalfluxes.hh  eulerfluxes.hh
mhd_eqns.hh  mhdfluxes.hh  mhd_eqns.cc  mhd_fluxes.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${fluxes_HEADERS} DESTINATION ${fluxesdir})