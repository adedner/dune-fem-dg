dune_install( models.hh problemcreator.hh steppertraits.hh)
     
set(MAININCLUDE_PATH  "${dune-fem-dg_SOURCE_DIR}/dune/fem-dg/main/")

include_directories(${MAININCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/advdiff/")

configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 

add_executable( advdiffall "${MAININCLUDE_PATH}/main.cc"  
                            "${MAININCLUDE_PATH}/main_0.cc" 
                            "${MAININCLUDE_PATH}/main_1.cc"
                            "${MAININCLUDE_PATH}/main_2.cc" 
                            "${MAININCLUDE_PATH}/main_3.cc")

dune_target_link_libraries(advdiffall "${DUNE_LIBS}")

add_executable( advdiff "${MAININCLUDE_PATH}/main.cc" "${MAININCLUDE_PATH}/main_pol.cc")
dune_target_link_libraries(advdiff "${DUNE_LIBS}")
set_property(TARGET advdiff APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=2")

set( PROG "advdiff" )

add_executable( ${PROG}_generatecode "${MAININCLUDE_PATH}/main.cc"  
                                     "${MAININCLUDE_PATH}/main_0.cc" 
                                     "${MAININCLUDE_PATH}/main_1.cc"
                                     "${MAININCLUDE_PATH}/main_2.cc" 
                                     "${MAININCLUDE_PATH}/main_3.cc")
set_property(TARGET ${PROG}_generatecode APPEND PROPERTY 
             COMPILE_DEFINITIONS "NDEBUG;BASEFUNCTIONSET_CODEGEN_GENERATE")
dune_target_link_libraries( ${PROG}_generatecode "${DUNE_LIBS}")

add_executable( ${PROG}_compilecode "${MAININCLUDE_PATH}/main.cc"  
                                     "${MAININCLUDE_PATH}/main_0.cc" 
                                     "${MAININCLUDE_PATH}/main_1.cc"
                                     "${MAININCLUDE_PATH}/main_2.cc" 
                                     "${MAININCLUDE_PATH}/main_3.cc")
set_property(TARGET ${PROG}_compilecode APPEND PROPERTY 
             COMPILE_DEFINITIONS "USE_BASEFUNCTIONSET_CODEGEN")
dune_target_link_libraries( ${PROG}_compilecode "${DUNE_LIBS}" )
           
           
#add_custom_command( OUTPUT generate_command
  #  COMMAND ${CMAKE_COMMAND} -D RUN_GENERATE_CODE_PROGRAM=${CMAKE_CURRENT_BINARY_DIR}/${PROG}_generatecode -D RUN_GENERATE_CODE_PARAMFILE="" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/RunGenerateCode.cmake )

  #  add_custom_command( OUTPUT generate_code_command
  #  COMMAND ${PROG}_generatecode
  #  COMMAND generate_command )
  #add_custom_target( generate_code DEPENDS generate_code_command )
  #
  #
  #add_custom_command( OUTPUT compile_code_command
  #  COMMAND ${PROG}_compilecode )
  #add_custom_target( compile_code DEPENDS compile_code_command )
  #
  #
  #add_custom_command( OUTPUT codegen_command
  #  COMMAND generate_code_command
  #  COMMAND compile_code_command )
  #add_custom_target( codegen DEPENDS codegen_command )



add_custom_target( generate 
  ${CMAKE_COMMAND} -D RUN_GENERATE_CODE_PROGRAM=${CMAKE_CURRENT_BINARY_DIR}/${PROG}_generatecode -D RUN_GENERATE_CODE_PARAMFILE="" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/RunGenerate.cmake ) 

add_custom_target( generate_code 
  ${CMAKE_COMMAND} -D RUN_GENERATE_CODE_EXEC=${PROG}_generatecode -D RUN_GENERATE_CODE_PROGRAM=${CMAKE_CURRENT_BINARY_DIR}/${PROG}_generatecode -D RUN_GENERATE_CODE_PARAMFILE="" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/RunGenerateCode.cmake ) 

add_custom_target( compile_code 
                   ${PROG}_compilecode )


add_custom_target( codegen
                   generate_code
                   COMMAND compile_code )
                 add_dependencies( codegen generate_code compile_code )

#add_dependencies( advdiff generate )

#generatecode:
#	$(MAKE) -i clean
#	$(MAKE) CODEGEN_OBJFILE=  CXXFLAGS="-g -Wall -Wfatal-errors -DBASEFUNCTIONSET_CODEGEN_GENERATE" $(PROG)
#	$(MAKE) generate
#
#generate:
#	./$(PROG) femhowto.eocSteps:1 femdg.stepper.maximaltimesteps:1 fem.io.outputformat:none fem.ode.order:1 paramFiles/paramNSWaves 
#
#compilecode:
#	$(MAKE) clean 
#	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DUSE_BASEFUNCTIONSET_CODEGEN" $(PROG)
#
#codegen:
#	$(MAKE) generatecode
#	$(MAKE) compilecode
#
#
dune_add_test( advdiff advdiffall )
