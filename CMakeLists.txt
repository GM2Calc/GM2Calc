cmake_minimum_required(VERSION 3.1)

project(GM2Calc VERSION 1.5.2)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED YES)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${PROJECT_SOURCE_DIR}/cmake/Eigen
  ${PROJECT_SOURCE_DIR}/cmake/Mathematica)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EXEC_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(GM2CALC_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include/gm2calc")
set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(INCLUDE_INSTALL_SUB_DIR "${INCLUDE_INSTALL_DIR}/gm2calc")
set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")

# generate gm2calc.pc file for pkg-config
configure_file (
  "${PROJECT_SOURCE_DIR}/cmake/gm2calc.pc.in"
  "${CMAKE_BINARY_DIR}/gm2calc.pc"
  @ONLY)

if(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup()
endif()

find_package(Boost 1.37.0 REQUIRED)
find_package(Eigen3 3.1 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(Mathematica COMPONENTS MathLink)

include_directories(${GM2CALC_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(doc)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
  enable_testing()
  add_subdirectory(test)
endif()

# make tag/release-* targets
find_package(Git)

if(Git_FOUND)
  add_custom_target(tag
    ${GIT_EXECUTABLE} tag v${PROJECT_VERSION}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Creating git tag v${PROJECT_VERSION}")

  add_custom_target(release-tag
    ${GIT_EXECUTABLE} archive --worktree-attributes --output=${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz --prefix=${PROJECT_NAME}-${PROJECT_VERSION}/ v${PROJECT_VERSION}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Creating ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz")

  add_custom_target(release-head
    ${GIT_EXECUTABLE} archive --worktree-attributes --output=${PROJECT_NAME}-HEAD.tar.gz --prefix=${PROJECT_NAME}-HEAD/ HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Creating ${PROJECT_NAME}-HEAD.tar.gz")
endif()

# install public headers
install(
    DIRECTORY ${GM2CALC_INCLUDE_DIRS}
    DESTINATION "${INCLUDE_INSTALL_DIR}"
  )

# install gm2calc.pc file
install(FILES
  "${CMAKE_BINARY_DIR}/gm2calc.pc"
  DESTINATION "${LIB_INSTALL_DIR}/pkgconfig")
