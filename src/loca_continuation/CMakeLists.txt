set(OBJS_LOCA_CONTINUATION
    loca_nln_problem.cpp loca_nln_group.cpp loca_nln_stepper.cpp loca_nln_statustest_factory.cpp
    )

baci_add_library(locacontinuation ${OBJS_LOCA_CONTINUATION})
target_link_libraries(locacontinuation PUBLIC nlnsolver_nox)
