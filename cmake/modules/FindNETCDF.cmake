mark_as_advanced(NETCDF_ROOT NETCDF_SUFFIX NETCDF_INCLUDEDIR NETCDF_LIBDIR)


find_path(NETCDF_INCLUDE_DIR
  NAMES "netcdf.h"
  PATHS ${NETCDF_ROOT} ${NETCDF_INCLUDEDIR}
  DOC "netcdf include directories")

# check header usability
include(CMakePushCheckState)
cmake_push_check_state()
set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} ${MPI_DUNE_COMPILE_FLAGS} -DENABLE_NETCDF")
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_DUNE_INCLUDE_PATH} ${NETCDF_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_DUNE_LIBRARIES})
include(CheckIncludeFiles)
check_include_files(netcdf.h NETCDF_HEADER_USABLE)

find_library(NETCDF_LIBRARY
  NAMES netcdf
  PATHS ${NETCDF_ROOT} ${NETCDF_LIBDIR}
  PATH_SUFFIXES "lib" "lib32" "lib64"
  DOC "netcdf library")

# check if library sion/sionser works
cmake_pop_check_state()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "NETCDF"
  REQUIRED_VARS NETCDF_INCLUDE_DIR NETCDF_LIBRARY
)

if (NETCDF_FOUND)
  set(NETCDF_INCLUDE_DIRS "${NETCDF_INCLUDE_DIR}")
  set(NETCDF_LIBRARIES "${NETCDF_LIBRARY}")

  #set HAVE_SIONLIB for config.h
  set(HAVE_NETCDF ${NETCDF_FOUND})

  dune_register_package_flags(LIBRARIES "${NETCDF_LIBRARIES}"
                              INCLUDE_DIRS "${NETCDF_INCLUDE_DIRS}")

endif ()
