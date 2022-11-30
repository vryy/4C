set(OBJS_DRT_VOLMORTAR
    volmortar_cell.cpp volmortar_coupling.cpp volmortar_integrator.cpp volmortar_utils.cpp
    )

baci_add_library(drt_volmortar ${OBJS_DRT_VOLMORTAR})
target_link_libraries(drt_volmortar PUBLIC drt_lib drt_lib_mortar)
