set(OBJS_DRT_CONSTRAINT
    constraint.cpp
    constraint_manager.cpp
    constraintdofset.cpp
    constraintpenalty.cpp
    constraintsolver.cpp
    monitor.cpp
    multipointconstraint.cpp
    multipointconstraint2.cpp
    multipointconstraint3.cpp
    multipointconstraint3penalty.cpp
    lagpenconstraint_noxinterface.cpp
    springdashpot_manager.cpp
    springdashpot.cpp
    nox_nln_lagpenconstraint_linearsystem.cpp
    )

baci_add_library(drt_constraint ${OBJS_DRT_CONSTRAINT})
target_link_libraries(drt_constraint PUBLIC linalg)

set(OBJS_DRT_LIB_CONSTRAINT
    constraint_element2.cpp
    constraint_element2_evaluate.cpp
    constraint_element3.cpp
    constraint_element3_evaluate.cpp
    )

baci_add_library(drt_lib_constraint ${OBJS_DRT_LIB_CONSTRAINT})
target_link_libraries(drt_lib_constraint PUBLIC drt_lib)
