set(OBJS_DRT_S8_LIB
    shell8.cpp
    shell8_evaluate.cpp
    shell8_input.cpp
    shell8_line.cpp
    shell8_line_evaluate.cpp
    shell8_scatra.cpp
    shell8/s8_mat_linel.c
    shell8/s8_mat_ogden.c
    shell8/s8_mat_ogden2.c
    shell8/s8_mat_ogden3.c
    shell8/s8_mattransform.c
    shell8/s8_tmat.c
    )

baci_add_library(ele_s8 ${OBJS_DRT_S8_LIB})
target_link_libraries(ele_s8 PUBLIC linalg)

baci_add_library(ele_s8_fortran fortran/s8_fserv.f)
target_link_libraries(ele_s8_fortran PRIVATE fortran)

target_link_libraries(ele_s8 PRIVATE ele_s8_fortran)
