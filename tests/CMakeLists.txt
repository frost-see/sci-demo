add_executable(test_derivative test_derivative.cpp)

target_link_libraries(test_derivative PRIVATE sci_derivative)

add_test(NAME test_derivative COMMAND test_derivative)
