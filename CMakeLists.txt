cmake_minimum_required(VERSION 3.15)
project(sci-demo VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)

add_subdirectory(src)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
