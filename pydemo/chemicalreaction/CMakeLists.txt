# add custom target to build tool chain for euler
dune_python_add_test(NAME chemical_python
                     SCRIPT ${PYTHON_EXECUTABLE} test.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
