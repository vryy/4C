set(OBJS_ELE_LUBRICATION
    # cmake-format: sortable
    lubrication_ele.cpp
    lubrication_ele_calc.cpp
    lubrication_ele_evaluate.cpp
    lubrication_ele_factory.cpp
    lubrication_ele_parameter.cpp
    )

baci_add_library(ele_lubrication ${OBJS_ELE_LUBRICATION})

set(LIBS_ELE_LUBRICATION
    # cmake-format: sortable
    lib
    )

baci_add_dependency(ele_lubrication ${LIBS_ELE_LUBRICATION})
