set(OBJS_W1_LIB
    # cmake-format: sortable
    wall1.cpp
    wall1_eas.cpp
    wall1_evaluate.cpp
    wall1_gemm.cpp
    wall1_input.cpp
    wall1_line.cpp
    wall1_line_evaluate.cpp
    wall1_mat.cpp
    wall1_nurbs.cpp
    wall1_poro.cpp
    wall1_poro_eletypes.cpp
    wall1_poro_evaluate.cpp
    wall1_poro_p1.cpp
    wall1_poro_p1_eletypes.cpp
    wall1_poro_p1_evaluate.cpp
    wall1_poro_p1_scatra.cpp
    wall1_poro_p1_scatra_eletypes.cpp
    wall1_poro_scatra.cpp
    wall1_poro_scatra_eletypes.cpp
    wall1_scatra.cpp
    wall1_scatra_evaluate.cpp
    )

baci_add_library(ele_w1 ${OBJS_W1_LIB})

set(LIBS_ELEW1
    # cmake-format: sortable
    lib
    lib_contact
    )

target_link_libraries(ele_w1 PUBLIC ${LIBS_ELEW1})
