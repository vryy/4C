set(OBJS_ELE_LUBRICATION
    lubrication_ele.cpp
    lubrication_ele_evaluate.cpp
    lubrication_ele_factory.cpp
    lubrication_ele_calc.cpp
    lubrication_ele_parameter.cpp
    )

baci_add_library(ele_lubrication ${OBJS_ELE_LUBRICATION})
target_link_libraries(ele_lubrication PUBLIC drt_lib)
