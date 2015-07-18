# the main files for all the tests are in the same dir
set(FEMDG_MAIN_DIR "${CMAKE_SOURCE_DIR}/dune/fem-dg/main")

# helper variables
set(SOURCEMAIN  ${FEMDG_MAIN_DIR}/main.cc)
set(SOURCEONE   ${SOURCEMAIN} ${FEMDG_MAIN_DIR}/main_pol.cc)
set(SOURCEALL   ${SOURCEMAIN} ${FEMDG_MAIN_DIR}/main_0.cc ${FEMDG_MAIN_DIR}/main_1.cc ${FEMDG_MAIN_DIR}/main_2.cc ${FEMDG_MAIN_DIR}/main_3.cc ${FEMDG_MAIN_DIR}/main_4.cc)

include(TargetDistclean)
message(AUTHOR_WARNING "TODO: Implement module test.")
