# Module that checks whether BGQ_L1PREFETCH is available.
#
# Variables used by this module which you may want to set:
# BGQ_L1PREFETCH_ROOT        Path list to search for BGQ_L1PREFETCH
#
# Sets the following variables
#
# BGQ_L1PREFETCH_FOUND          True if BGQ_L1PREFETCH was found and usable
# HAVE_BGQ_L1PREFETCH           True if BGQ_L1PREFETCH was found and usable
# BGQ_L1PREFETCH_INCLUDE_DIRS   Path to the BGQ_L1PREFETCH include dirs
# BGQ_L1PREFETCH_LIBRARIES      Name of the BGQ_L1PREFETCH libraries
#

set(BGQ_L1PREFETCH_ROOT "" CACHE PATH "Path list to search for BGQ_L1PREFETCH")


#look for header files at positions given by the user
find_path(BGQ_L1PREFETCH_INCLUDE_DIR
  NAME "l1p/pprefetch.h"
  PATHS ${BGQ_L1PREFETCH_ROOT}
  PATH_SUFFIXES "bgq_l1prefetch"
  NO_DEFAULT_PATH
)
#now also look for default paths
find_path(BGQ_L1PREFETCH_INCLUDE_DIR
  NAMES "l1p/pprefetch.h"
  PATH_SUFFIXES "bgq_l1prefetch" 
)

set(BGQ_L1PREFETCH_INCLUDEDIR "${BGQ_L1PREFETCH_INCLUDE_DIR}/include" CACHE PATH "directory with BGQ_L1PREFETCH headers inside")
set(BGQ_L1PREFETCH_LIBDIR "${BGQ_L1PREFETCH_INCLUDE_DIR}/lib" CACHE PATH "directory with BGQ_L1PREFETCH libraries inside")

mark_as_advanced(BGQ_L1PREFETCH_ROOT BGQ_L1PREFETCH_INCLUDEDIR BGQ_L1PREFETCH_LIBDIR)

# check header usability
include(CMakePushCheckState)
cmake_push_check_state()
set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -DENABLE_BGQ_L1PREFETCH=1")
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${BGQ_L1PREFETCH_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES} ${PARMETIS_LIBRARIES} -L${BGQ_L1PREFETCH_LIBDIR}")
check_include_files(zoltan_cpp.h BGQ_L1PREFETCH_HEADER_USABLE)

#look for library at positions given by the user
find_library(BGQ_L1PREFETCH_LIBRARY
  NAMES "SPI_l1p"
  PATHS ${BGQ_L1PREFETCH_LIBDIR}
  NO_DEFAULT_PATH
)

# check if library zoltan works
include(CheckSymbolExists)
if(BGQ_L1PREFETCH_LIBRARY)
  get_filename_component(BGQ_L1PREFETCH_LIB_PATH ${BGQ_L1PREFETCH_LIBRARY} PATH)
  check_library_exists(SPI_lp1 main ${BGQ_L1PREFETCH_LIB_PATH} BGQ_L1PREFETCH_LIB_WORKS)
endif(BGQ_L1PREFETCH_LIBRARY)

cmake_pop_check_state()


# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "BGQ_L1PREFETCH"
  DEFAULT_MSG
  BGQ_L1PREFETCH_INCLUDE_DIR
  BGQ_L1PREFETCH_LIBRARY
  BGQ_L1PREFETCH_HEADER_USABLE
  BGQ_L1PREFETCH_LIB_WORKS
)

mark_as_advanced(BGQ_L1PREFETCH_INCLUDE_DIR BGQ_L1PREFETCH_LIBRARY BGQ_L1PREFETCH_LIB_WORKS BGQ_L1PREFETCH_HEADER_USABLE)

# if both headers and library are found, store results
if(BGQ_L1PREFETCH_FOUND)
  set(BGQ_L1PREFETCH_INCLUDE_DIRS ${BGQ_L1PREFETCH_INCLUDE_DIR})
  set(BGQ_L1PREFETCH_LIBRARIES ${BGQ_L1PREFETCH_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of BGQ_L1PREFETCH succeded:\n"
    "Include directory: ${BGQ_L1PREFETCH_INCLUDE_DIRS}\n"
    "Library directory: ${BGQ_L1PREFETCH_LIBRARIES}\n\n")
  set(BGQ_L1PREFETCH_DUNE_COMPILE_FLAGS "-I${BGQ_L1PREFETCH_INCLUDE_DIRS} -DENABLE_BGQ_L1PREFETCH=1"
    CACHE STRING "Compile Flags used by DUNE when compiling with BGQ_L1PREFETCH programs")
  set(BGQ_L1PREFETCH_DUNE_LIBRARIES ${BGQ_L1PREFETCH_LIBRARIES} 
    CACHE STRING "Libraries used by DUNE when linking BGQ_L1PREFETCH programs")
else(BGQ_L1PREFETCH_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of BGQ_L1PREFETCH failed:\n"
    "Include directory: ${BGQ_L1PREFETCH_INCLUDE_DIRS}\n"
    "Library directory: ${BGQ_L1PREFETCH_LIBRARIES}\n\n")
endif(BGQ_L1PREFETCH_FOUND)

#set HAVE_BGQ_L1PREFETCH for config.h
set(HAVE_BGQ_L1PREFETCH ${BGQ_L1PREFETCH_FOUND})

#add all bgq_l1prefetch related flags to ALL_PKG_FLAGS, this must happen regardless of a target using add_dune_bgq_l1prefetch_flags
if(BGQ_L1PREFETCH_FOUND)
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "${BGQ_L1PREFETCH_DUNE_COMPILE_FLAGS}")
  foreach(dir "${BGQ_L1PREFETCH_INCLUDE_DIRS}")
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
  endforeach()
endif()
