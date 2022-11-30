set(OBJS_DRT_SURFSTRESS drt_surfstress_manager.cpp)

baci_add_library(drt_lib_surfstress ${OBJS_DRT_SURFSTRESS})
target_link_libraries(drt_lib_surfstress PUBLIC drt_lib)
