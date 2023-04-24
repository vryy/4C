set(OBJS_LOCA_CONTINUATION
    # cmake-format: sortable
    loca_continuation_nln_group.cpp
    loca_continuation_nln_problem.cpp
    loca_continuation_nln_statustest_factory.cpp
    loca_continuation_nln_stepper.cpp
    )

baci_add_library(locacontinuation ${OBJS_LOCA_CONTINUATION})

set(LIBS_LOCACONTINUATION
    # cmake-format: sortable
    nlnsolver_nox
    )

baci_add_dependency(locacontinuation ${LIBS_LOCACONTINUATION})
