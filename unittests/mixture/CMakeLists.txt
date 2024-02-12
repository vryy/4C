set(TESTNAME unittests_mixture)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_mixture_full_constrained_mixture_fiber_adaptive_history_test.cpp
    baci_mixture_full_constrained_mixture_fiber_test.cpp
    baci_remodelfiber_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_target_requires_cyclic_legacy_linking(${TESTNAME})
