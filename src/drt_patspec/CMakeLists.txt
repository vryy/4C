set(OBJS_DRT_PATSPEC patspec.cpp)

baci_add_library(drt_patspec ${OBJS_DRT_PATSPEC})
target_link_libraries(drt_patspec PUBLIC drt_mat)
