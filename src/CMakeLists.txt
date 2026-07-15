add_library(sci_derivative STATIC derivative.cpp)
add_library(sci_integration STATIC integration.cpp)

target_include_directories(sci_derivative PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(sci_integration PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
