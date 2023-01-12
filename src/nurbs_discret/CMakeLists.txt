set(OBJS_NURBS_DISCRET
    # cmake-format: sortable
    nurbs_discret.cpp
    nurbs_discret_apply_nurbs_initial_condition.cpp
    nurbs_discret_control_point.cpp
    nurbs_discret_knotvector.cpp
    )

baci_add_library(nurbs_discret ${OBJS_NURBS_DISCRET})

set(LIBS_NURBS_DISCRET
    # cmake-format: sortable
    lib
    )

target_link_libraries(nurbs_discret PUBLIC ${LIBS_NURBS_DISCRET})
