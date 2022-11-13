set(OBJS_SOLVER
    bgs2x2_operator.cpp
    solver_cheapsimplepreconditioner.cpp
    solver_directsolver.cpp
    solver_krylovsolver.cpp
    solver_pointpreconditioner.cpp
    solver_blockpreconditioners.cpp
    solver_krylovprojectionpreconditioner.cpp
    solver_ifpackpreconditioner.cpp
    solver_mlpreconditioner.cpp
    solver_muelupreconditioner.cpp
    solver_belossolver.cpp
    solver_aztecsolver.cpp
    solver_aztecsolver_projectedresidual.cpp
    solver_amgnxn_preconditioner.cpp
    solver_amgnxn_vcycle.cpp
    solver_amgnxn_hierarchies.cpp
    solver_amgnxn_smoothers.cpp
    solver_amgnxn_objects.cpp
    linalg_solver.cpp
    )

baci_add_library(solver ${OBJS_SOLVER})
target_link_libraries(solver PUBLIC linalg)
