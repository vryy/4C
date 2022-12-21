set(OBJS_VOLMORTAR
    # cmake-format: sortable
    volmortar_cell.cpp
    volmortar_coupling.cpp
    volmortar_integrator.cpp
    volmortar_utils.cpp
    )

baci_add_library(volmortar ${OBJS_VOLMORTAR})

set(LIBS_VOLMORTAR
    # cmake-format: sortable
    lib
    lib_mortar
    mortar
    )

target_link_libraries(volmortar PUBLIC ${LIBS_VOLMORTAR})
