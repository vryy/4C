set(OBJS_DRT_TO3_LIB
    # cmake-format: sortable
    torsion3.cpp
    torsion3_evaluate.cpp
    torsion3_input.cpp
    )

baci_add_library(ele_to3 ${OBJS_DRT_TO3_LIB})

set(LIBS_ELE_TO3
    # cmake-format: sortable
    drt_lib
    drt_lib_browniandyn
    )

target_link_libraries(ele_to3 PUBLIC ${LIBS_ELE_TO3})
