# add custom target to build tool chain for euler
#add_custom_target(euler_python
#                  COMMAND ${PYTHON_EXECUTABLE} main.py
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# add custom target to build tool chain for euler
add_test(NAME euler_python
         COMMAND ${PYTHON_EXECUTABLE} test.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
