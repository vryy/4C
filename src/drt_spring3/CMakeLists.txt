set(OBJS_DRT_SP3_LIB spring3.cpp spring3_evaluate.cpp spring3_input.cpp)

baci_add_library(ele_sp3 ${OBJS_DRT_SP3_LIB})
target_link_libraries(ele_sp3 PUBLIC drt_lib)
