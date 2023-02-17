set(OBJS_TUTORIAL
    # cmake-format: sortable
    tutorial_dyn.cpp
    tutorial_fixedpoint_scheme.cpp
    tutorial_inpar.cpp
    tutorial_nln_truss.cpp
    tutorial_struct_ale_coupling.cpp
    tutorial_utils.cpp
    )

baci_add_library(alg_tutorial ${OBJS_TUTORIAL})

set(LIBS_ALG_TUTORIAL
    # cmake-format: sortable
    alg_tutorial_ele
    alg_tutorial_material
    inpar
    linalg
    )

baci_add_dependency(alg_tutorial ${LIBS_ALG_TUTORIAL})

add_subdirectory(tutorial_ele)
add_subdirectory(tutorial_material)
