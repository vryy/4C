set(OBJS_FIBER_LIB
    # cmake-format: sortable
    fiber_nodal_fiber_holder.cpp
    fiber_node.cpp
    fiber_utils.cpp
    )

baci_add_library(lib_fiber ${OBJS_FIBER_LIB})

set(LIBS_LIB_FIBER
    # cmake-format: sortable
    lib
    linalg
    )

baci_add_dependency(lib_fiber ${LIBS_LIB_FIBER})
