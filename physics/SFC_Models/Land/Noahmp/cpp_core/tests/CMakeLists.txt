# GFS Noah-MP C++ Core test suite
cmake_minimum_required(VERSION 3.18)

add_executable(test_noahmpdrv_parity test_noahmpdrv_parity.F90)
target_include_directories(test_noahmpdrv_parity PRIVATE ..)
target_link_libraries(test_noahmpdrv_parity PRIVATE ccpp_physics)
add_test(NAME test_noahmpdrv_parity COMMAND test_noahmpdrv_parity)

add_executable(benchmark_noahmpdrv benchmark_noahmpdrv.cpp)
target_include_directories(benchmark_noahmpdrv PRIVATE ..)
target_link_libraries(benchmark_noahmpdrv PRIVATE ccpp_physics)
add_test(NAME benchmark_noahmpdrv COMMAND benchmark_noahmpdrv)
