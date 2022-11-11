set(OBJS_DRT_SOLID_ELE
    solid_ele_calc.cpp
    solid_ele_evaluate.cpp
    solid_ele_factory.cpp
    solid_ele.cpp
    solid_utils.cpp
    )

baci_add_library(drt_solid_ele ${OBJS_DRT_SOLID_ELE})
baci_link_default_libraries(drt_solid_ele)

# TODO: This module requires more includes that are currently included transitively from drt_lib
target_link_libraries(drt_solid_ele PUBLIC drt_lib)
