# GM2Calc library
add_library(gm2calc
  gm2_1loop_c.cpp
  gm2_1loop.cpp
  gm2_2loop_c.cpp
  gm2_2loop.cpp
  gm2_dilog.cpp
  gm2_ffunctions.cpp
  gm2_mb.cpp
  gm2_numerics.cpp
  gm2_slha_io.cpp
  gm2_uncertainty_c.cpp
  gm2_uncertainty.cpp
  MSSMNoFV_onshell_c.cpp
  MSSMNoFV_onshell.cpp
  MSSMNoFV_onshell_mass_eigenstates.cpp
  MSSMNoFV_onshell_physical.cpp
  MSSMNoFV_onshell_problems.cpp
  MSSMNoFV_onshell_soft_parameters.cpp
  MSSMNoFV_onshell_susy_parameters.cpp
  )
target_include_directories(gm2calc
  PUBLIC
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_BINARY_DIR}/include"
    "${EIGEN3_INCLUDE_DIR}"
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_library(GM2Calc::GM2Calc ALIAS gm2calc)

install(
  TARGETS gm2calc
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# GM2Calc main executable
add_executable(gm2calc.x gm2calc.cpp)
target_link_libraries(gm2calc.x PRIVATE GM2Calc::GM2Calc)
target_include_directories(gm2calc.x PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

install(
  TARGETS gm2calc.x
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# MathML executable
if(Mathematica_FOUND)
  Mathematica_MathLink_ADD_EXECUTABLE(
    gm2calc.mx
    gm2calc.tm)
  target_link_libraries(gm2calc.mx GM2Calc::GM2Calc)
  target_include_directories(gm2calc.mx
    PRIVATE ${Mathematica_MathLink_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  install(
    TARGETS gm2calc.mx
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
