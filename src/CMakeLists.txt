# GM2Calc library
add_library(gm2calc
  gm2_dilog.cpp
  gm2_ffunctions.cpp
  gm2_mb.cpp
  gm2_numerics.cpp
  gm2_slha_io.cpp
  MSSMNoFV/gm2_1loop_c.cpp
  MSSMNoFV/gm2_1loop.cpp
  MSSMNoFV/gm2_2loop_c.cpp
  MSSMNoFV/gm2_2loop.cpp
  MSSMNoFV/gm2_uncertainty_c.cpp
  MSSMNoFV/gm2_uncertainty.cpp
  MSSMNoFV/MSSMNoFV_onshell_c.cpp
  MSSMNoFV/MSSMNoFV_onshell.cpp
  MSSMNoFV/MSSMNoFV_onshell_mass_eigenstates.cpp
  MSSMNoFV/MSSMNoFV_onshell_physical.cpp
  MSSMNoFV/MSSMNoFV_onshell_problems.cpp
  MSSMNoFV/MSSMNoFV_onshell_soft_parameters.cpp
  MSSMNoFV/MSSMNoFV_onshell_susy_parameters.cpp
  GeneralTHDM/gm2_1loop.cpp
  GeneralTHDM/gm2_1loop_H.cpp
  GeneralTHDM/gm2_2loop.cpp
  GeneralTHDM/gm2_2loop_B.cpp
  GeneralTHDM/gm2_2loop_F.cpp
  GeneralTHDM/GeneralTHDM.cpp
  GeneralTHDM/GeneralTHDM_mass_eigenstates.cpp
  GeneralTHDM/GeneralTHDM_physical.cpp
  GeneralTHDM/GeneralTHDM_problems.cpp
  GeneralTHDM/GeneralTHDM_soft_parameters.cpp
  GeneralTHDM/GeneralTHDM_susy_parameters.cpp
  )
target_include_directories(gm2calc
  PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    "${EIGEN3_INCLUDE_DIR}"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${Boost_INCLUDE_DIRS}"
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
if(Mathematica_MathLink_FOUND)
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
