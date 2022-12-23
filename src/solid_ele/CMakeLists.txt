set(OBJS_ELE_SOLID
    solid_ele_calc.cpp
    solid_ele_evaluate.cpp
    solid_ele_factory.cpp
    solid_ele.cpp
    solid_utils.cpp
    )

baci_add_library(ele_solid ${OBJS_ELE_SOLID})
baci_link_default_libraries(ele_solid)

# TODO: This module requires more includes that are currently included transitively from lib
target_link_libraries(ele_solid PUBLIC lib)
