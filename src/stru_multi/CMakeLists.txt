set(OBJS_STRU_MULTI
    # cmake-format: sortable
    stru_multi_convert_mat.cpp
    stru_multi_microstatic.cpp
    stru_multi_microstatic_npsupport.cpp
    stru_multi_microstatic_service.cpp
    )

baci_add_library(stru_multi ${OBJS_STRU_MULTI})

set(LIBS_STRU_MULTI
    # cmake-format: sortable
    inpar
    lib
    )

target_link_libraries(stru_multi PUBLIC ${LIBS_STRU_MULTI})
