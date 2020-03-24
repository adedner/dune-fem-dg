# add custom target to build tool chain for euler
dune_python_add_test(NAME euler_pythonA
         COMMAND ${PYTHON_EXECUTABLE} test.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
dune_python_add_test(NAME euler_pythonB
         COMMAND ${PYTHON_EXECUTABLE} test.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
dune_python_add_test(NAME euler_pythonC
         COMMAND ${PYTHON_EXECUTABLE} test.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
dune_python_add_test(NAME euler_pythonD
         COMMAND ${PYTHON_EXECUTABLE} test.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
dune_python_add_test(NAME euler_pythonE
         COMMAND ${PYTHON_EXECUTABLE} test.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
dune_python_add_test(NAME euler_pythonF
         COMMAND ${PYTHON_EXECUTABLE} test.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
dune_python_add_test(NAME euler_python
         COMMAND ${PYTHON_EXECUTABLE} test.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
