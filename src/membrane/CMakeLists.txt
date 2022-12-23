set(OBJS_MEMBRANE
    # cmake-format: sortable
    membrane.cpp
    membrane_eletypes.cpp
    membrane_evaluate.cpp
    membrane_input.cpp
    membrane_line.cpp
    membrane_line_evaluate.cpp
    membrane_scatra.cpp
    membrane_scatra_eletypes.cpp
    membrane_scatra_evaluate.cpp
    )

baci_add_library(ele_membrane ${OBJS_MEMBRANE})

set(LIBS_ELE_MEMBRANE
    # cmake-format: sortable
    lib
    )

target_link_libraries(ele_membrane PUBLIC ${LIBS_ELE_MEMBRANE})
