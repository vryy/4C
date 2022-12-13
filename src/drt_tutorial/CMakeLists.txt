set(OBJS_DRT_TUTORIAL
    # cmake-format: sortable
    inpar_tutorial.cpp
    tutorial_dyn.cpp
    tutorial_fixedpoint_scheme.cpp
    tutorial_nln_truss.cpp
    tutorial_struct_ale_coupling.cpp
    tutorial_utils.cpp
    )

baci_add_library(alg_tutorial ${OBJS_DRT_TUTORIAL})

set(LIBS_ALG_TUTORIAL
    # cmake-format: sortable
    alg_tutorial_ele
    alg_tutorial_material
    drt_inpar
    linalg
    )

target_link_libraries(alg_tutorial PUBLIC ${LIBS_ALG_TUTORIAL})

add_subdirectory(tutorial_ele)
add_subdirectory(tutorial_material)
