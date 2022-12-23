set(OBJS_TUTORIAL_MATERIAL
    # cmake-format: sortable
    tutorial_mat_stvenantkirchhoff.cpp
    )

baci_add_library(alg_tutorial_material ${OBJS_TUTORIAL_MATERIAL})

set(LIBS_ALG_TUTORIAL_MATERIAL
    # cmake-format: sortable
    alg_tutorial
    )

target_link_libraries(alg_tutorial_material PUBLIC ${LIBS_ALG_TUTORIAL_MATERIAL})
