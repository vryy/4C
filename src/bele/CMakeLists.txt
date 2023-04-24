set(OBJS_ELE_BELE
    # cmake-format: sortable
    bele_bele2.cpp
    bele_bele2_evaluate.cpp
    bele_bele3.cpp
    bele_bele3_evaluate.cpp
    bele_bele3_line.cpp
    bele_bele3_line_evaluate.cpp
    bele_vele3.cpp
    bele_vele3_evaluate.cpp
    bele_vele3_line.cpp
    bele_vele3_line_evaluate.cpp
    bele_vele3_surface.cpp
    bele_vele3_surface_evaluate.cpp
    )

baci_add_library(ele_bele ${OBJS_ELE_BELE})

set(LIBS_ELE_BELE
    # cmake-format: sortable
    comm
    ele_scatra
    ele_so
    lib
    )

baci_add_dependency(ele_bele ${LIBS_ELE_BELE})
