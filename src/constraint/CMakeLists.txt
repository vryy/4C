set(OBJS_CONSTRAINT
    # cmake-format: sortable
    constraint.cpp
    constraint_dofset.cpp
    constraint_lagpenconstraint_noxinterface.cpp
    constraint_manager.cpp
    constraint_monitor.cpp
    constraint_multipointconstraint.cpp
    constraint_multipointconstraint2.cpp
    constraint_multipointconstraint3.cpp
    constraint_multipointconstraint3penalty.cpp
    constraint_nox_nln_lagpenconstraint_linearsystem.cpp
    constraint_penalty.cpp
    constraint_solver.cpp
    constraint_springdashpot.cpp
    constraint_springdashpot_manager.cpp
    )

baci_add_library(constraint ${OBJS_CONSTRAINT})

set(LIBS_CONSTRAINT
    # cmake-format: sortable
    linalg
    )

baci_add_dependency(constraint ${LIBS_CONSTRAINT})

set(OBJS_LIB_CONSTRAINT
    # cmake-format: sortable
    constraint_element2.cpp
    constraint_element2_evaluate.cpp
    constraint_element3.cpp
    constraint_element3_evaluate.cpp
    )

baci_add_library(lib_constraint ${OBJS_LIB_CONSTRAINT})

set(LIBS_LIB_CONSTRAINT
    # cmake-format: sortable
    lib
    )

baci_add_dependency(lib_constraint ${LIBS_LIB_CONSTRAINT})
