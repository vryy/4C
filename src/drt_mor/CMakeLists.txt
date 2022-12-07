set(OBJS_DRT_MOR mor_pod.cpp)

baci_add_library(drt_mor ${OBJS_DRT_MOR})
target_link_libraries(drt_mor PUBLIC linalg)
