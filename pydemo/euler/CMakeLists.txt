dune_fem_add_elliptic_models(euler.ufl)

add_executable( pyeuler EXCLUDE_FROM_ALL main.cc euler.hh )
target_compile_definitions( pyeuler PRIVATE YASPGRID;GRIDDIM=2;POLORDER=2 )
add_dune_mpi_flags( pyeuler )
target_link_libraries( pyeuler )
