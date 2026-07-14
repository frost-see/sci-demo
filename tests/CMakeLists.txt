include(FetchContent)
include(GoogleTest)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)

add_executable(test_derivative test_derivative.cpp)

target_link_libraries(test_derivative PRIVATE sci_derivative GTest::gtest_main)

gtest_discover_tests(test_derivative)

add_executable(test_curve_fitting test_curve_fitting.cpp)

target_link_libraries(test_curve_fitting PRIVATE sci_curve_fitting GTest::gtest_main)

gtest_discover_tests(test_curve_fitting)
