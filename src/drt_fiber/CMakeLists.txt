set(OBJS_DRT_FIBER_LIB
    # cmake-format: sortable
    drt_fiber_node.cpp
    drt_fiber_utils.cpp
    nodal_fiber_holder.cpp
    )

baci_add_library(drt_lib_fiber ${OBJS_DRT_FIBER_LIB})

set(LIBS_DRT_LIB_FIBER
    # cmake-format: sortable
    drt_lib
    linalg
    )

target_link_libraries(drt_lib_fiber PUBLIC ${LIBS_DRT_LIB_FIBER})
