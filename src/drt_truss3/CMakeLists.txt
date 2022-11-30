set(OBJS_DRT_T3_LIB truss3.cpp truss3_evaluate.cpp truss3_input.cpp truss3_scatra.cpp)

baci_add_library(ele_t3 ${OBJS_DRT_T3_LIB})
target_link_libraries(ele_t3 PUBLIC drt_lib)
