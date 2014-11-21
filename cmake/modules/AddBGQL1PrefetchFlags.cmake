# Module providing convenience functions for using 
#
# Provides the following functions:
#
# add_dune_bgq_l1prefetch_flags(target1 target2 ...)
#
# Adds the necessary flags to compile and link the targets with BGQ_L1PREFETCH support.
#
function(add_dune_bgq_l1prefetch_flags _targets)
  if(BGQ_L1PREFETCH_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} ${BGQ_L1PREFETCH_LIBRARIES})
    endforeach(_target ${_targets})
    set_property(TARGET ${_targets}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS ENABLE_BGQ_L1PREFETCH=1 )
    set_property(TARGET ${_targets} APPEND PROPERTY
      INCLUDE_DIRECTORIES "${BGQ_L1PREFETCH_INCLUDE_DIRS}")
  endif(BGQ_L1PREFETCH_FOUND)
endfunction(add_dune_bdq_l1prefetch_flags _targets)
