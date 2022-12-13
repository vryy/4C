set(OBJS_DRT_SP3_LIB
    # cmake-format: sortable
    spring3.cpp
    spring3_evaluate.cpp
    spring3_input.cpp
    )

baci_add_library(ele_sp3 ${OBJS_DRT_SP3_LIB})

set(LIBS_ELE_SP3
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(ele_sp3 PUBLIC ${LIBS_ELE_SP3})
