set(OBJS_TUTORIAL_MATERIAL
    # cmake-format: sortable
    tutorial_material_stvenantkirchhoff.cpp
    )

baci_add_library(alg_tutorial_material ${OBJS_TUTORIAL_MATERIAL})

set(LIBS_ALG_TUTORIAL_MATERIAL
    # cmake-format: sortable
    alg_tutorial
    )

baci_add_dependency(alg_tutorial_material ${LIBS_ALG_TUTORIAL_MATERIAL})
