# add custom target distclean
#if(${CMAKE_GENERATOR} MATCHES ".*Unix Makefiles.*")
#add_custom_target(distclean @echo cleaning CMake files)
#set(DISTCLEANFILES
#cmake.depends
#cmake.check_depends
#CMakeCache.txt
#*/CMakeCache.txt
#CMakeFiles
#*/CMakeFiles
#Makefile
#*/Makefile
#)

#add_custom_command(
#DEPENDS clean
#COMMENT "distribution clean"
#COMMAND rm
#ARGS -Rf ${DISTCLEANFILES}
#TARGET distclean
#)
#endif()
