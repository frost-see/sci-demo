include(FetchContent)
include(GoogleTest)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)

add_executable(test_derivative test_derivative.cpp)
add_executable(test_integration test_integration.cpp)

target_link_libraries(test_derivative PRIVATE sci_derivative GTest::gtest_main)
target_link_libraries(test_integration PRIVATE sci_integration GTest::gtest_main)

gtest_discover_tests(test_derivative)
gtest_discover_tests(test_integration)
