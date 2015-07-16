
message( "${RUN_GENERATE_CODE_PROGRAM} fem.eoc.steps:1 femdg.stepper.maximaltimesteps:1 fem.io.outputformat:none ${RUN_GENERATE_CODE_PARAMFILE} " )
execute_process(COMMAND 
  ${RUN_GENERATE_CODE_PROGRAM} fem.eoc.steps:1 femdg.stepper.maximaltimesteps:1 fem.io.outputformat:none ${RUN_GENERATE_CODE_PARAMFILE} )

