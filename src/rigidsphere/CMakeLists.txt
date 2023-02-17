set(OBJS_ELE_RIGIDSPHERE
    # cmake-format: sortable
    rigidsphere.cpp
    rigidsphere_evaluate.cpp
    rigidsphere_input.cpp
    )

baci_add_library(ele_rigidsphere ${OBJS_ELE_RIGIDSPHERE})

set(LIBS_ELE_RIGIDSPHERE
    # cmake-format: sortable
    linalg
    )

baci_add_dependency(ele_rigidsphere ${LIBS_ELE_RIGIDSPHERE})
