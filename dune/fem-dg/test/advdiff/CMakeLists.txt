dune_install( models.hh problemcreator.hh steppertraits.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/advdiff/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 

add_executable( advdiffall main.cc )
dune_target_link_libraries(advdiffall "${DUNE_LIBS}")

add_executable( advdiff main.cc )
dune_target_link_libraries(advdiff "${DUNE_LIBS}")
set_property(TARGET advdiff APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=2")

           #add_code_generate_targets( advdiff main.cc )
           #add_code_generate_targets( advdiffall main.cc )
           #
           #add_custom_target( generate 
           #  ${CMAKE_COMMAND} -D RUN_GENERATE_CODE_PROGRAM=${CMAKE_CURRENT_BINARY_DIR}/${PROG}_generatecode -D RUN_GENERATE_CODE_PARAMFILE="" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/RunGenerate.cmake ) 
           #
           #add_custom_target( generate_code 
           #  ${CMAKE_COMMAND} -D RUN_GENERATE_CODE_EXEC=${PROG}_generatecode -D RUN_GENERATE_CODE_PROGRAM=${CMAKE_CURRENT_BINARY_DIR}/${PROG}_generatecode -D RUN_GENERATE_CODE_PARAMFILE="" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/RunGenerateCode.cmake ) 
           #
           #add_custom_target( compile_code 
           #                   ${PROG}_compilecode )
           #
           #
           #add_custom_target( codegen
           #                   generate_code
           #                   COMMAND compile_code )
           #
           #add_dependencies( codegen generate_code compile_code )

#add_dependencies( advdiff generate )

#generatecode:
#	$(MAKE) -i clean
#	$(MAKE) CODEGEN_OBJFILE=  CXXFLAGS="-g -Wall -Wfatal-errors -DBASEFUNCTIONSET_CODEGEN_GENERATE" $(PROG)
#	$(MAKE) generate
#
#generate:
#	./$(PROG) fem.eoc.steps:1 femdg.stepper.maximaltimesteps:1 fem.io.outputformat:none fem.ode.order:1 paramFiles/paramNSWaves 
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
