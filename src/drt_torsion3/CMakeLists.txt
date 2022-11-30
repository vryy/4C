set(OBJS_DRT_TO3_LIB torsion3.cpp torsion3_evaluate.cpp torsion3_input.cpp)

baci_add_library(ele_to3 ${OBJS_DRT_TO3_LIB})
target_link_libraries(ele_to3 PUBLIC drt_lib drt_lib_browniandyn)
