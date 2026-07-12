add_library(sci_derivative STATIC derivative.cpp)

target_include_directories(sci_derivative PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(sci_curve_fitting STATIC curve_fitting.cpp)

target_include_directories(sci_curve_fitting PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
