set(OBJS_DRT_NURBS_DISCRET
    drt_apply_nurbs_initial_condition.cpp
    drt_control_point.cpp
    drt_knotvector.cpp
    drt_nurbs_discret.cpp
    )

baci_add_library(drt_nurbs_discret ${OBJS_DRT_NURBS_DISCRET})
target_link_libraries(drt_nurbs_discret PUBLIC drt_lib)
