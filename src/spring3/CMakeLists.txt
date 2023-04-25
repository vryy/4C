set(OBJS_SP3_LIB
    # cmake-format: sortable
    spring3.cpp
    spring3_evaluate.cpp
    spring3_input.cpp
    )

baci_add_library(ele_sp3 ${OBJS_SP3_LIB})

set(LIBS_ELE_SP3
    # cmake-format: sortable
    lib
    )

baci_add_dependency(ele_sp3 ${LIBS_ELE_SP3})
