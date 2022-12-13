set(OBJS_DRT_NURBS_DISCRET
    # cmake-format: sortable
    drt_apply_nurbs_initial_condition.cpp
    drt_control_point.cpp
    drt_knotvector.cpp
    drt_nurbs_discret.cpp
    )

baci_add_library(drt_nurbs_discret ${OBJS_DRT_NURBS_DISCRET})

set(LIBS_DRT_NURBS_DISCRET
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(drt_nurbs_discret PUBLIC ${LIBS_DRT_NURBS_DISCRET})
