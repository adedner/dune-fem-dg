dune_install(
  main.cc
  algorithmcreator.hh
  models.hh
  navierstokesflux.hh
  problems.hh
  thermodynamics.hh
)
dune_add_subdirs( test problems )
