set(OBJS_SOLVER_NONLIN_NOX
    # cmake-format: sortable
    solver_nonlin_nox_aux.cpp
    solver_nonlin_nox_constraint_group.cpp
    solver_nonlin_nox_direction_defaultsteptest.cpp
    solver_nonlin_nox_direction_factory.cpp
    solver_nonlin_nox_direction_modified_newton.cpp
    solver_nonlin_nox_direction_newton.cpp
    solver_nonlin_nox_floating_point_exception.cpp
    solver_nonlin_nox_globaldata.cpp
    solver_nonlin_nox_group.cpp
    solver_nonlin_nox_group_prepostoperator.cpp
    solver_nonlin_nox_inner_statustest_armijo.cpp
    solver_nonlin_nox_inner_statustest_combo.cpp
    solver_nonlin_nox_inner_statustest_factory.cpp
    solver_nonlin_nox_inner_statustest_filter.cpp
    solver_nonlin_nox_inner_statustest_generic.cpp
    solver_nonlin_nox_inner_statustest_upperbound.cpp
    solver_nonlin_nox_inner_statustest_volume_change.cpp
    solver_nonlin_nox_linearsystem.cpp
    solver_nonlin_nox_linearsystem_factory.cpp
    solver_nonlin_nox_linearsystem_prepostoperator.cpp
    solver_nonlin_nox_linesearch_backtrack.cpp
    solver_nonlin_nox_linesearch_factory.cpp
    solver_nonlin_nox_linesearch_prepostoperator.cpp
    solver_nonlin_nox_meritfunction_factory.cpp
    solver_nonlin_nox_meritfunction_infeasibility.cpp
    solver_nonlin_nox_meritfunction_lagrangian.cpp
    solver_nonlin_nox_problem.cpp
    solver_nonlin_nox_solver_factory.cpp
    solver_nonlin_nox_solver_linesearchbased.cpp
    solver_nonlin_nox_solver_prepostop_generic.cpp
    solver_nonlin_nox_solver_ptc.cpp
    solver_nonlin_nox_statustest_activeset.cpp
    solver_nonlin_nox_statustest_combo.cpp
    solver_nonlin_nox_statustest_factory.cpp
    solver_nonlin_nox_statustest_normf.cpp
    solver_nonlin_nox_statustest_normupdate.cpp
    solver_nonlin_nox_statustest_normwrms.cpp
    )

baci_add_library(nlnsolver_nox ${OBJS_SOLVER_NONLIN_NOX})

set(LIBS_NLNSOLVER_NOX
    # cmake-format: sortable
    linalg
    )

target_link_libraries(nlnsolver_nox PUBLIC ${LIBS_NLNSOLVER_NOX})
