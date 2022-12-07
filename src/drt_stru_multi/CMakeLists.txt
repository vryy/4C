set(OBJS_DRT_STRU_MULTI
    convert_mat.cpp microstatic.cpp microstatic_npsupport.cpp microstatic_service.cpp
    )

baci_add_library(drt_stru_multi ${OBJS_DRT_STRU_MULTI})
target_link_libraries(drt_stru_multi PUBLIC drt_lib drt_inpar)
