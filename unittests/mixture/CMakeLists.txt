set(TESTNAME unittests_mixture)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_mixture_full_constrained_mixture_fiber_adaptive_history_test.cpp
    4C_mixture_full_constrained_mixture_fiber_test.cpp
    4C_remodelfiber_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
