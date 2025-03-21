add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

function(add_gm2calc_test name ext)
  add_executable(${name}.x ${name}.${ext})
  target_link_libraries(${name}.x PRIVATE GM2Calc::GM2Calc)
  target_include_directories(${name}.x PRIVATE $<TARGET_PROPERTY:GM2Calc::GM2Calc,INCLUDE_DIRECTORIES>)
  if(MSVC)
    target_compile_options(${name}.x PRIVATE /bigobj)
  endif()
  add_test(NAME ${name} COMMAND ${name}.x)
endfunction()

function(add_gm2calc_bench name ext)
  add_executable(${name}.x ${name}.${ext})
  target_link_libraries(${name}.x PRIVATE GM2Calc::GM2Calc)
  target_include_directories(${name}.x PRIVATE $<TARGET_PROPERTY:GM2Calc::GM2Calc,INCLUDE_DIRECTORIES>)
endfunction()

add_gm2calc_bench(test_benchmark           cpp)
add_gm2calc_bench(test_benchmark_ffunctions cpp)
add_gm2calc_test(test_dilog                cpp)
add_gm2calc_test(test_eigen_utils          cpp)
add_gm2calc_test(test_ffunctions           cpp)
add_gm2calc_test(test_linalg_1             cpp)
add_gm2calc_test(test_linalg_2             cpp)
add_gm2calc_test(test_linalg_3             cpp)
add_gm2calc_test(test_linalg_4             cpp)
add_gm2calc_test(test_linalg_5             cpp)
add_gm2calc_test(test_linalg_6             cpp)
add_gm2calc_test(test_linalg_7             cpp)
add_gm2calc_test(test_linalg_8             cpp)
add_gm2calc_test(test_mf                   cpp)
add_gm2calc_test(test_MSSMNoFV             cpp)
add_gm2calc_test(test_MSSMNoFV_c_interface cpp)
add_gm2calc_test(test_MSSMNoFV_slha_io     cpp)
add_gm2calc_test(test_numerics             cpp)
add_gm2calc_test(test_SM                   cpp)
add_gm2calc_test(test_SM_c_interface       cpp)
add_gm2calc_test(test_SM_slha_io           cpp)
add_gm2calc_test(test_slha_io              cpp)
add_gm2calc_test(test_THDM                 cpp)
add_gm2calc_test(test_THDM_amu             cpp)
add_gm2calc_test(test_THDM_c_interface     cpp)
add_gm2calc_test(test_THDM_slha_io         cpp)
add_gm2calc_test(test_version              cpp)

# test shell scripts
find_program(BASH_PROGRAM sh)

if(BASH_PROGRAM)
  add_test(NAME test_c_examples
    COMMAND ${BASH_PROGRAM} -c
    "BINDIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${CMAKE_CURRENT_SOURCE_DIR}/test_c_examples.sh > test_c_examples.sh.out 2>&1")
  add_test(NAME test_points
    COMMAND ${BASH_PROGRAM} -c
    "GM2CALC=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gm2calc.x ${CMAKE_CURRENT_SOURCE_DIR}/test_points.sh > test_points.sh.out 2>&1")
  add_test(NAME test_gm2calc
    COMMAND ${BASH_PROGRAM} -c
    "GM2CALC=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gm2calc.x ${CMAKE_CURRENT_SOURCE_DIR}/test_gm2calc.sh > test_gm2calc.sh.out 2>&1")
endif (BASH_PROGRAM)

# test Mathematica scripts
if(Mathematica_MathLink_FOUND)
  Mathematica_MathLink_ADD_TEST(
    NAME test_mathlink_interface
    TARGET gm2calc.mx
    SCRIPT test_mathlink_interface.m)
endif()

# test Python scripts
if(Python_INTERFACE)
  configure_file(${PROJECT_SOURCE_DIR}/test/test_THDM_python_interface.py.in 
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_THDM_python_interface.py 
                 @ONLY)
  add_test(
    NAME test_THDM_python_interface
    COMMAND "${PYTHON_EXE}" test_THDM_python_interface.py
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  configure_file(${PROJECT_SOURCE_DIR}/test/test_MSSMNoFV_python_interface.py.in 
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_MSSMNoFV_python_interface.py 
                 @ONLY)
  add_test(
    NAME test_MSSMNoFV_python_interface
    COMMAND "${PYTHON_EXE}" test_MSSMNoFV_python_interface.py
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()
