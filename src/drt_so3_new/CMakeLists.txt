set(OBJS_DRT_SOLID_ELE_NEW
    solid_ele_calc.cpp
    solid_ele_evaluate.cpp
    solid_ele_factory.cpp
    solid_ele.cpp
    solid_utils.cpp
    )

baci_add_library(drt_solid_ele_new ${OBJS_DRT_SOLID_ELE_NEW})
baci_link_default_libraries(drt_solid_ele_new)

# TODO: This module requires more includes that are currently included transitively from drt_lib
target_link_libraries(drt_solid_ele_new PUBLIC drt_lib)
