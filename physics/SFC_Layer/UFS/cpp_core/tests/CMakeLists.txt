# GFS Surface Layer C++ Core test suite
cmake_minimum_required(VERSION 3.18)

add_executable(test_sfc_diff_parity test_sfc_diff_parity.F90)
target_include_directories(test_sfc_diff_parity PRIVATE ..)
target_link_libraries(test_sfc_diff_parity PRIVATE ccpp_physics)
add_test(NAME test_sfc_diff_parity COMMAND test_sfc_diff_parity)

add_executable(test_sfc_nst_parity test_sfc_nst_parity.F90)
target_include_directories(test_sfc_nst_parity PRIVATE ..)
target_link_libraries(test_sfc_nst_parity PRIVATE ccpp_physics)
add_test(NAME test_sfc_nst_parity COMMAND test_sfc_nst_parity)

add_executable(benchmark_sfc_physics benchmark_sfc_physics.cpp)
target_include_directories(benchmark_sfc_physics PRIVATE ..)
target_link_libraries(benchmark_sfc_physics PRIVATE ccpp_physics)
add_test(NAME benchmark_sfc_physics COMMAND benchmark_sfc_physics)
