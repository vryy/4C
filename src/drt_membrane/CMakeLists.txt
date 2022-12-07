set(OBJS_DRT_MEMBRANE
    membrane_eletypes.cpp
    membrane.cpp
    membrane_evaluate.cpp
    membrane_input.cpp
    membrane_line.cpp
    membrane_line_evaluate.cpp
    membrane_scatra_eletypes.cpp
    membrane_scatra.cpp
    membrane_scatra_evaluate.cpp
    )

baci_add_library(ele_membrane ${OBJS_DRT_MEMBRANE})
target_link_libraries(ele_membrane PUBLIC drt_lib)
