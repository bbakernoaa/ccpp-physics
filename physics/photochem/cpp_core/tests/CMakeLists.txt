# GFS Photochemistry C++ Core test suite
cmake_minimum_required(VERSION 3.18)

add_executable(test_ozphys_parity test_ozphys_parity.F90)
target_include_directories(test_ozphys_parity PRIVATE ..)
target_link_libraries(test_ozphys_parity PRIVATE ccpp_physics)
add_test(NAME test_ozphys_parity COMMAND test_ozphys_parity)

add_executable(test_h2ophys_parity test_h2ophys_parity.F90)
target_include_directories(test_h2ophys_parity PRIVATE ..)
target_link_libraries(test_h2ophys_parity PRIVATE ccpp_physics)
add_test(NAME test_h2ophys_parity COMMAND test_h2ophys_parity)

add_executable(benchmark_photochem benchmark_photochem.cpp)
target_include_directories(benchmark_photochem PRIVATE ..)
target_link_libraries(benchmark_photochem PRIVATE ccpp_physics)
add_test(NAME benchmark_photochem COMMAND benchmark_photochem)
