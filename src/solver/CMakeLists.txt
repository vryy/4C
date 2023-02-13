set(OBJS_SOLVER
    # cmake-format: sortable
    solver_amgnxn_hierarchies.cpp
    solver_amgnxn_objects.cpp
    solver_amgnxn_preconditioner.cpp
    solver_amgnxn_smoothers.cpp
    solver_amgnxn_vcycle.cpp
    solver_belossolver.cpp
    solver_bgs2x2_operator.cpp
    solver_blockpreconditioners.cpp
    solver_cheapsimplepreconditioner.cpp
    solver_directsolver.cpp
    solver_ifpackpreconditioner.cpp
    solver_krylovprojectionpreconditioner.cpp
    solver_krylovsolver.cpp
    solver_linalg_ana.cpp
    solver_linalg_precond.cpp
    solver_linalg_solver.cpp
    solver_mlpreconditioner.cpp
    solver_muelupreconditioner.cpp
    solver_pointpreconditioner.cpp
    )

baci_add_library(solver ${OBJS_SOLVER})

set(LIBS_SOLVER
    # cmake-format: sortable
    linalg
    )

target_link_libraries(solver PUBLIC ${LIBS_SOLVER})
