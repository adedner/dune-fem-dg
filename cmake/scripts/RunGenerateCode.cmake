
message("COMMAND ${RUN_GENERATE_CODE_EXEC}" )
execute_process(COMMAND make ${RUN_GENERATE_CODE_EXEC} )
message( "${RUN_GENERATE_CODE_PROGRAM} fem.eoc.steps:1 femdg.stepper.maximaltimesteps:1 fem.io.outputformat:none fem.ode.order:1 ${RUN_GENERATE_CODE_PARAMFILE} " )
execute_process(COMMAND 
  ${RUN_GENERATE_CODE_PROGRAM} fem.eoc.steps:1 femdg.stepper.maximaltimesteps:1 fem.io.outputformat:none fem.ode.order:1 ${RUN_GENERATE_CODE_PARAMFILE} )

