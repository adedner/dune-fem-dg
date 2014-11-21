set(stepperdir  ${top_srcdir}/dune/fem-dg/stepper/)
set(stepper_HEADERS  advectiondiffusionstepper.hh
                  advectionstepper.hh
                  baseevolution.hh
                  base.hh
                  ellipt.hh
                  stepperbase.hh
                  stepper.hh
                  steppertraits.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${stepper_HEADERS} DESTINATION ${stepperdir})