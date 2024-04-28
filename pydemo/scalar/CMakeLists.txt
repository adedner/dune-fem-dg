# add custom target to build tool chain for scalar problems
dune_python_add_test(NAME scalardg_python
                     SCRIPT testdg.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(scalardg_python PROPERTIES TIMEOUT 2500)
