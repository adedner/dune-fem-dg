dune_python_add_test(NAME bndtest_python
                     SCRIPT bndtest.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(bndtest_python PROPERTIES TIMEOUT 2500)
