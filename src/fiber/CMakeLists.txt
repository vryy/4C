set(OBJS_FIBER_LIB
    # cmake-format: sortable
    fiber_node.cpp
    fiber_utils.cpp
    nodal_fiber_holder.cpp
    )

baci_add_library(lib_fiber ${OBJS_FIBER_LIB})

set(LIBS_LIB_FIBER
    # cmake-format: sortable
    lib
    linalg
    )

target_link_libraries(lib_fiber PUBLIC ${LIBS_LIB_FIBER})
