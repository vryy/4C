set(OBJS_W1_LIB
    # cmake-format: sortable
    w1.cpp
    w1_eas.cpp
    w1_evaluate.cpp
    w1_gemm.cpp
    w1_input.cpp
    w1_line.cpp
    w1_line_evaluate.cpp
    w1_mat.cpp
    w1_nurbs.cpp
    w1_poro.cpp
    w1_poro_eletypes.cpp
    w1_poro_evaluate.cpp
    w1_poro_p1.cpp
    w1_poro_p1_eletypes.cpp
    w1_poro_p1_evaluate.cpp
    w1_poro_p1_scatra.cpp
    w1_poro_p1_scatra_eletypes.cpp
    w1_poro_scatra.cpp
    w1_poro_scatra_eletypes.cpp
    w1_scatra.cpp
    w1_scatra_evaluate.cpp
    )

baci_add_library(ele_w1 ${OBJS_W1_LIB})

set(LIBS_ELEW1
    # cmake-format: sortable
    lib
    lib_contact
    )

baci_add_dependency(ele_w1 ${LIBS_ELEW1})
