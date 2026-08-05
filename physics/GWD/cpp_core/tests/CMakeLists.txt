# GWD C++ Core test suite
cmake_minimum_required(VERSION 3.18)

add_executable(test_gwd_parity test_gwd_parity.F90)
target_include_directories(test_gwd_parity PRIVATE ..)
target_link_libraries(test_gwd_parity PRIVATE ccpp_physics)
add_test(NAME test_gwd_parity COMMAND test_gwd_parity)

add_executable(benchmark_gwd benchmark_gwd.cpp)
target_include_directories(benchmark_gwd PRIVATE ..)
target_link_libraries(benchmark_gwd PRIVATE ccpp_physics)
add_test(NAME benchmark_gwd COMMAND benchmark_gwd)
