configure_file(config.h.in config.h)

add_library(gm2calc
  dilog.cpp
  ffunctions.cpp
  gm2_1loop_c.cpp
  gm2_1loop.cpp
  gm2_2loop_c.cpp
  gm2_2loop.cpp
  gm2_mb.cpp
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
target_include_directories(gm2calc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# GM2Calc main executable
add_executable(gm2calc.x gm2calc.cpp)
target_link_libraries(gm2calc.x gm2calc)
target_include_directories(gm2calc.x PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# MathML executable
