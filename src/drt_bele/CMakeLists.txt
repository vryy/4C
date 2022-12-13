set(OBJS_ELE_BELE
    # cmake-format: sortable
    bele2.cpp
    bele2_evaluate.cpp
    bele3.cpp
    bele3_evaluate.cpp
    bele3_line.cpp
    bele3_line_evaluate.cpp
    vele3.cpp
    vele3_evaluate.cpp
    vele3_line.cpp
    vele3_line_evaluate.cpp
    vele3_surface.cpp
    vele3_surface_evaluate.cpp
    )

baci_add_library(ele_bele ${OBJS_ELE_BELE})

set(LIBS_ELE_BELE
    # cmake-format: sortable
    drt_comm
    drt_lib
    ele_scatra
    ele_so
    )

target_link_libraries(ele_bele PUBLIC ${LIBS_ELE_BELE})
