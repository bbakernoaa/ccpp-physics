# SATMEDMF C++ Core test suite

add_executable(test_implicit_stability test_implicit_stability.cpp)
target_include_directories(test_implicit_stability PRIVATE ..)
target_link_libraries(test_implicit_stability PRIVATE ccpp_physics)
add_test(NAME test_implicit_stability COMMAND test_implicit_stability)

add_executable(test_pbl_lifecycle test_pbl_lifecycle.cpp)
target_include_directories(test_pbl_lifecycle PRIVATE ..)
target_link_libraries(test_pbl_lifecycle PRIVATE ccpp_physics)
add_test(NAME test_pbl_lifecycle COMMAND test_pbl_lifecycle)

add_executable(test_pbl_strides test_pbl_strides.cpp)
target_include_directories(test_pbl_strides PRIVATE ..)
target_link_libraries(test_pbl_strides PRIVATE ccpp_physics)
add_test(NAME test_pbl_strides COMMAND test_pbl_strides)

add_executable(test_local_mixing test_local_mixing.cpp)
target_include_directories(test_local_mixing PRIVATE ..)
target_link_libraries(test_local_mixing PRIVATE ccpp_physics)
add_test(NAME test_local_mixing COMMAND test_local_mixing)

add_executable(test_updrafts test_updrafts.cpp)
target_include_directories(test_updrafts PRIVATE ..)
target_link_libraries(test_updrafts PRIVATE ccpp_physics)
add_test(NAME test_updrafts COMMAND test_updrafts)

add_executable(test_downdrafts test_downdrafts.cpp)
target_include_directories(test_downdrafts PRIVATE ..)
target_link_libraries(test_downdrafts PRIVATE ccpp_physics)
add_test(NAME test_downdrafts COMMAND test_downdrafts)

add_executable(test_canopy test_canopy.cpp)
target_include_directories(test_canopy PRIVATE ..)
target_link_libraries(test_canopy PRIVATE ccpp_physics)
add_test(NAME test_canopy COMMAND test_canopy)

add_executable(test_canopy_mass_conservation test_canopy_mass_conservation.cpp)
target_include_directories(test_canopy_mass_conservation PRIVATE ..)
target_link_libraries(test_canopy_mass_conservation PRIVATE ccpp_physics)
add_test(NAME test_canopy_mass_conservation COMMAND test_canopy_mass_conservation)

add_executable(test_satmedmf_parity test_satmedmf_parity.F90)
target_include_directories(test_satmedmf_parity PRIVATE ..)
target_link_libraries(test_satmedmf_parity PRIVATE ccpp_physics)
add_test(NAME test_satmedmf_parity COMMAND test_satmedmf_parity)

add_executable(benchmark_satmedmf_cpp benchmark_satmedmf_cpp.cpp)
target_include_directories(benchmark_satmedmf_cpp PRIVATE ..)
target_link_libraries(benchmark_satmedmf_cpp PRIVATE ccpp_physics)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../satmedmfvdifq.F ${CMAKE_CURRENT_BINARY_DIR}/satmedmfvdifq_fortran.F COPYONLY)

add_executable(benchmark_satmedmf
               benchmark_satmedmf.F90
               ${CMAKE_CURRENT_BINARY_DIR}/satmedmfvdifq_fortran.F)
target_compile_definitions(benchmark_satmedmf PRIVATE "TESTING_PARITY")
target_include_directories(benchmark_satmedmf PRIVATE ..)
target_link_libraries(benchmark_satmedmf PRIVATE ccpp_physics)

