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

baci_add_dependency(nurbs_discret ${LIBS_NURBS_DISCRET})
