set(OBJS_DRT_TUTORIAL
    inpar_tutorial.cpp
    tutorial_dyn.cpp
    tutorial_fixedpoint_scheme.cpp
    tutorial_nln_truss.cpp
    tutorial_struct_ale_coupling.cpp
    tutorial_utils.cpp
    )

add_subdirectory(tutorial_ele)
add_subdirectory(tutorial_material)

baci_add_library(alg_tutorial ${OBJS_DRT_TUTORIAL})
target_link_libraries(
  alg_tutorial PUBLIC alg_tutorial_material alg_tutorial_ele drt_lib_inpar linalg
  )
