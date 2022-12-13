set(OBJS_DRT_VOLMORTAR
    # cmake-format: sortable
    volmortar_cell.cpp
    volmortar_coupling.cpp
    volmortar_integrator.cpp
    volmortar_utils.cpp
    )

baci_add_library(drt_volmortar ${OBJS_DRT_VOLMORTAR})

set(LIBS_DRT_VOLMORTAR
    # cmake-format: sortable
    drt_lib
    drt_lib_mortar
    drt_mortar
    )

target_link_libraries(drt_volmortar PUBLIC ${LIBS_DRT_VOLMORTAR})
