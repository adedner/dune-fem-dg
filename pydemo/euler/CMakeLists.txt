# add custom target to build tool chain for euler
dune_python_add_test(NAME eulerdg_python
                     SCRIPT testdg.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(eulerdg_python PROPERTIES TIMEOUT 2500)

dune_python_add_test(NAME eulerfv1_python
                     SCRIPT testfv1.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(eulerfv1_python PROPERTIES TIMEOUT 2500)

dune_python_add_test(NAME eulerfv0_python
                     SCRIPT testfv0.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(eulerfv0_python PROPERTIES TIMEOUT 2500)
