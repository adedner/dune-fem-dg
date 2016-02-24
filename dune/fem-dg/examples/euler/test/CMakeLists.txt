include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/euler/")

# configure parameter file from parameter.in
configure_parameter_file()

# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )

if( NOT FEMDG_FAST_TESTBUILD )
  if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
    add_executable( eulerall ../main.cc ../problems/chorjo.cc )
    dune_target_enable_all_packages( eulerall )
    dune_add_test( eulerall )
  else()
    dune_add_test( NAME eulerall
                   SOURCES ../main.cc ../problems/chorjo.cc
                   COMPILE_DEFINITIONS LIMITER )
  endif()
  add_code_generate_targets( eulerall )
endif()

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  add_executable( euler ../main.cc ../problems/chorjo.cc )
  set_property(TARGET euler APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )
  dune_target_enable_all_packages( euler )
  dune_add_test( euler )
else()
  dune_add_test( NAME euler
                 SOURCES ../main.cc ../problems/chorjo.cc
                 COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )
endif()
add_code_generate_targets( euler )
