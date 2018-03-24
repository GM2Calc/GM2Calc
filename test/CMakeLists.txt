# test C/C++ executables
function(add_gm2calc_test name ext)
  add_executable(${name}.x ${name}.${ext})
  target_link_libraries(${name}.x GM2Calc::gm2calc)
  add_test(NAME ${name} COMMAND ${name}.x)
endfunction()

add_gm2calc_test(test_c_interface cpp)

# test shell scripts
find_program(BASH_PROGRAM bash)

if(BASH_PROGRAM)
  add_test(NAME test_c_examples
    COMMAND ${BASH_PROGRAM} -c
    "BINDIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${CMAKE_CURRENT_SOURCE_DIR}/test_c_examples.sh > test_c_examples.sh.out 2>&1")
  add_test(NAME test_points
    COMMAND ${BASH_PROGRAM} -c
    "GM2CALC=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gm2calc.x ${CMAKE_CURRENT_SOURCE_DIR}/test_points.sh > test_points.sh.out 2>&1")
endif (BASH_PROGRAM)

# test Mathematica scripts
if(Mathematica_FOUND)
  Mathematica_MathLink_ADD_TEST(
    NAME test_mathlink_interface
    TARGET gm2calc.mx
    SCRIPT test_mathlink_interface.m)
endif()
