set(OBJS_NURBS_DISCRET
    # cmake-format: sortable
    apply_nurbs_initial_condition.cpp
    control_point.cpp
    knotvector.cpp
    nurbs_discret.cpp
    )

baci_add_library(nurbs_discret ${OBJS_NURBS_DISCRET})

set(LIBS_NURBS_DISCRET
    # cmake-format: sortable
    lib
    )

target_link_libraries(nurbs_discret PUBLIC ${LIBS_NURBS_DISCRET})
