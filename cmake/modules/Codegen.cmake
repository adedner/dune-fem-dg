function(add_code_generate_targets _target)

  get_target_property( compiledefs ${_target} COMPILE_DEFINITIONS )
  string(REPLACE "compiledefs-NOTFOUND" "" compiledefs "${compiledefs}")

  #Add targets and executables
  add_executable( ${_target}_generatecode ${ARGN} )
  set_property(TARGET ${_target}_generatecode APPEND PROPERTY COMPILE_DEFINITIONS "NDEBUG;BASEFUNCTIONSET_CODEGEN_GENERATE;${compiledefs}")
  dune_target_link_libraries( ${_target}_generatecode "${DUNE_LIBS}" )

  add_executable( ${_target}_compilecode ${ARGN} )
  set_property(TARGET ${_target}_compilecode APPEND PROPERTY COMPILE_DEFINITIONS "NUSE_BASEFUNCTIONSET_CODEGEN;${compiledefs}")
  dune_target_link_libraries( ${_target}_compilecode "${DUNE_LIBS}" )

  add_custom_target( ${_target}_generate 
    ${CMAKE_COMMAND} -D RUN_CODEGEN_PROGRAM=${CMAKE_CURRENT_BINARY_DIR}/${_target}_generatecode -D RUN_CODEGEN_PARAMFILE="" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/RunGenerate.cmake ) 
  
  add_custom_target( ${_target}_codegen )

  #Add depenencies  
  add_dependencies( ${_target}_generate ${_target}_generatecode )
  add_dependencies( ${_target}_compilecode ${_target}_generate )
  add_dependencies( ${_target}_codegen ${_target}_compilecode )

endfunction(add_code_generate_targets _target)
