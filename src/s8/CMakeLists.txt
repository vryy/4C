set(OBJS_S8_LIB
    # cmake-format: sortable
    s8.cpp
    s8_evaluate.cpp
    s8_input.cpp
    s8_line.cpp
    s8_line_evaluate.cpp
    s8_scatra.cpp
    shell8/s8_shell8_mat_linel.c
    shell8/s8_shell8_mattransform.c
    shell8/s8_shell8_tmat.c
    )

baci_add_library(ele_s8 ${OBJS_S8_LIB})

set(LIBS_ELE_S8
    # cmake-format: sortable
    fortran
    linalg
    )

baci_add_dependency(ele_s8 ${LIBS_ELE_S8})
