cmake_minimum_required(VERSION 3.1)

project(GM2Calc VERSION 1.4.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Mathematica)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(Mathematica COMPONENTS MathLink)

enable_testing()
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(test)
add_subdirectory(doc)
