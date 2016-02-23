include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/euler/")

# configure parameter file from parameter.in
configure_parameter_file()

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
#add_definitions( "-DUSE_PTHREADS" )

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
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_target_enable_all_packages( euler )
  dune_add_test( euler )
else()
  dune_add_test( NAME euler
                 SOURCES ../main.cc ../problems/chorjo.cc
                 COMPILE_DEFINITIONS POLORDER=2 )
endif()
add_code_generate_targets( euler )
