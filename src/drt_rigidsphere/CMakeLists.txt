set(OBJS_DRT_RIGIDSPHERE_LIB rigidsphere.cpp rigidsphere_evaluate.cpp rigidsphere_input.cpp)

baci_add_library(ele_rigidsphere ${OBJS_DRT_RIGIDSPHERE_LIB})
target_link_libraries(ele_rigidsphere PUBLIC linalg)
