function(add_gm2calc_test name ext)
  add_executable(${name}.x ${name}.${ext})
  target_link_libraries(${name}.x gm2calc)
  add_test(NAME ${name} COMMAND ${name}.x)
endfunction()

add_gm2calc_test(test_c_interface cpp)
