set(OBJS_LOCA_CONTINUATION
    # cmake-format: sortable
    loca_nln_group.cpp
    loca_nln_problem.cpp
    loca_nln_statustest_factory.cpp
    loca_nln_stepper.cpp
    )

baci_add_library(locacontinuation ${OBJS_LOCA_CONTINUATION})

set(LIBS_LOCACONTINUATION
    # cmake-format: sortable
    nlnsolver_nox
    )

target_link_libraries(locacontinuation PUBLIC ${LIBS_LOCACONTINUATION})
