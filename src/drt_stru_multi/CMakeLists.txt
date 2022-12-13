set(OBJS_DRT_STRU_MULTI
    # cmake-format: sortable
    convert_mat.cpp
    microstatic.cpp
    microstatic_npsupport.cpp
    microstatic_service.cpp
    )

baci_add_library(drt_stru_multi ${OBJS_DRT_STRU_MULTI})

set(LIBS_DRT_STRU_MULTI
    # cmake-format: sortable
    drt_inpar
    drt_lib
    )

target_link_libraries(drt_stru_multi PUBLIC ${LIBS_DRT_STRU_MULTI})
