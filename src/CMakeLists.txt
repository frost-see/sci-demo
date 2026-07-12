add_library(sci_derivative STATIC derivative.cpp)

target_include_directories(sci_derivative PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
