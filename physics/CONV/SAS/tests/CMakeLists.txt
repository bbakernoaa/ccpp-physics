enable_testing()

add_executable(test_tracer_transport_all
    test_tracer_transport_all.F90
    ../tracer_transport_mod.F90
    ../../../hooks/machine.F
    ../../../hooks/physcons.F90
)
add_test(NAME test_tracer_transport_all COMMAND test_tracer_transport_all)

add_executable(test_multi_tracer
    test_multi_tracer.F90
    ../tracer_transport_mod.F90
    ../../../hooks/machine.F
    ../../../hooks/physcons.F90
)
add_test(NAME test_multi_tracer COMMAND test_multi_tracer)

add_executable(test_sas_comprehensive
    test_sas_comprehensive.F90
    ../sascnvn.F
    ../shalcnv.F
    ../tracer_transport_mod.F90
    ../../../hooks/machine.F
    ../../../hooks/physcons.F90
    ../../../tools/funcphys.f90
)
add_test(NAME test_sas_comprehensive COMMAND test_sas_comprehensive)

add_executable(test_sas_conservation
    test_sas_conservation.F90
    ../sascnvn.F
    ../shalcnv.F
    ../tracer_transport_mod.F90
    ../../../hooks/machine.F
    ../../../hooks/physcons.F90
    ../../../tools/funcphys.f90
)
add_test(NAME test_sas_conservation COMMAND test_sas_conservation)
