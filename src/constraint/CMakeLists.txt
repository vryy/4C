set(OBJS_CONSTRAINT
    # cmake-format: sortable
    constraint.cpp
    constraint_manager.cpp
    constraintdofset.cpp
    constraintpenalty.cpp
    constraintsolver.cpp
    lagpenconstraint_noxinterface.cpp
    monitor.cpp
    multipointconstraint.cpp
    multipointconstraint2.cpp
    multipointconstraint3.cpp
    multipointconstraint3penalty.cpp
    nox_nln_lagpenconstraint_linearsystem.cpp
    springdashpot.cpp
    springdashpot_manager.cpp
    )

baci_add_library(constraint ${OBJS_CONSTRAINT})

set(LIBS_CONSTRAINT
    # cmake-format: sortable
    linalg
    )

target_link_libraries(constraint PUBLIC ${LIBS_CONSTRAINT})

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

target_link_libraries(lib_constraint PUBLIC ${LIBS_LIB_CONSTRAINT})
