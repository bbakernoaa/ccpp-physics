# RRTMG C++ Core test suite
cmake_minimum_required(VERSION 3.18)

add_executable(test_rrtmg_sw_parity test_rrtmg_sw_parity.F90)
target_include_directories(test_rrtmg_sw_parity PRIVATE ..)
target_link_libraries(test_rrtmg_sw_parity PRIVATE ccpp_physics)
add_test(NAME test_rrtmg_sw_parity COMMAND test_rrtmg_sw_parity)

add_executable(test_rrtmg_lw_parity test_rrtmg_lw_parity.F90)
target_include_directories(test_rrtmg_lw_parity PRIVATE ..)
target_link_libraries(test_rrtmg_lw_parity PRIVATE ccpp_physics)
add_test(NAME test_rrtmg_lw_parity COMMAND test_rrtmg_lw_parity)

add_executable(benchmark_rrtmg_radiation benchmark_rrtmg_radiation.cpp)
target_include_directories(benchmark_rrtmg_radiation PRIVATE ..)
target_link_libraries(benchmark_rrtmg_radiation PRIVATE ccpp_physics)
add_test(NAME benchmark_rrtmg_radiation COMMAND benchmark_rrtmg_radiation)
