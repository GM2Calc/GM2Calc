cmake_minimum_required(VERSION 3.1)

project(GM2Calc
  DESCRIPTION "Calculation of the anomalous magnetic moment of the muon in the MSSM"
  LANGUAGES C CXX)

include(GNUInstallDirs)

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
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if(EXISTS ${PROJECT_BINARY_DIR}/conanbuildinfo.cmake)
  include(${PROJECT_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup()
endif()

find_package(Boost 1.37.0 REQUIRED)
find_package(Eigen3 3.1 REQUIRED)
find_package(Mathematica COMPONENTS MathLink)

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(doc)

if(EXISTS "${PROJECT_SOURCE_DIR}/test")
  enable_testing()
  add_subdirectory(test)
endif()

# install public headers
install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/include/gm2calc"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

# generate gm2calc.pc file for pkg-config
configure_file (
  "${PROJECT_SOURCE_DIR}/cmake/gm2calc.pc.in"
  "${PROJECT_BINARY_DIR}/gm2calc.pc"
  @ONLY)

# install gm2calc.pc file
install(FILES
  "${PROJECT_BINARY_DIR}/gm2calc.pc"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig")
