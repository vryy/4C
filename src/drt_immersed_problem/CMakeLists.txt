set(OBJS_DRT_IMMERSED_LIB
    # cmake-format: sortable
    immersed_node.cpp
    )

baci_add_library(drt_lib_immersed ${OBJS_DRT_IMMERSED_LIB})

set(LIBS_DRT_LIB_IMMERSED
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(drt_lib_immersed PUBLIC ${LIBS_DRT_LIB_IMMERSED})

set(OBJS_ALG_IMMERSED_PROBLEMS
    # cmake-format: sortable
    fsi_partitioned_immersed.cpp
    immersed_base.cpp
    immersed_partitioned_fsi_dirichletneumann.cpp
    immersed_problem_dyn.cpp
    )

baci_add_library(alg_immersed_problem ${OBJS_ALG_IMMERSED_PROBLEMS})

set(LIBS_ALG_IMMERSED_PROBLEMS
    # cmake-format: sortable
    alg_fsi
    )

target_link_libraries(alg_immersed_problem PUBLIC ${LIBS_ALG_IMMERSED_PROBLEMS})
