cmake_minimum_required(VERSION 3.1)

project(GM2Calc VERSION 1.4.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

enable_testing()
add_subdirectory(src)
add_subdirectory(examples)
# add_subdirectory(test)
