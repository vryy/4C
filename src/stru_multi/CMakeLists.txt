set(OBJS_STRU_MULTI
    # cmake-format: sortable
    convert_mat.cpp
    microstatic.cpp
    microstatic_npsupport.cpp
    microstatic_service.cpp
    )

baci_add_library(stru_multi ${OBJS_STRU_MULTI})

set(LIBS_STRU_MULTI
    # cmake-format: sortable
    inpar
    lib
    )

target_link_libraries(stru_multi PUBLIC ${LIBS_STRU_MULTI})
