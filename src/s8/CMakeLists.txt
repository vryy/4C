set(OBJS_S8_LIB
    # cmake-format: sortable
    s8.cpp
    s8_evaluate.cpp
    s8_input.cpp
    s8_line.cpp
    s8_line_evaluate.cpp
    s8_scatra.cpp
    shell8/s8_shell8_mat_linel.c
    shell8/s8_shell8_mat_ogden.c
    shell8/s8_shell8_mat_ogden2.c
    shell8/s8_shell8_mat_ogden3.c
    shell8/s8_shell8_mattransform.c
    shell8/s8_shell8_tmat.c
    )

baci_add_library(ele_s8 ${OBJS_S8_LIB})

set(LIBS_ELE_S8
    # cmake-format: sortable
    fortran
    linalg
    )

target_link_libraries(ele_s8 PUBLIC ${LIBS_ELE_S8})

baci_add_library(ele_s8_fortran fortran/s8_fserv.f)
target_link_libraries(ele_s8_fortran PRIVATE fortran)

target_link_libraries(ele_s8 PRIVATE ele_s8_fortran)
