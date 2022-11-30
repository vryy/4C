set(OBJS_DRT_FIBER_LIB drt_fiber_node.cpp drt_fiber_utils.cpp nodal_fiber_holder.cpp)

baci_add_library(drt_lib_fiber ${OBJS_DRT_FIBER_LIB})
target_link_libraries(drt_lib_fiber PUBLIC linalg drt_lib)
