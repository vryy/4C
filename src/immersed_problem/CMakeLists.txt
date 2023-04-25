set(OBJS_IMMERSED_LIB
    # cmake-format: sortable
    immersed_problem_immersed_node.cpp
    )

baci_add_library(lib_immersed ${OBJS_IMMERSED_LIB})

set(LIBS_LIB_IMMERSED
    # cmake-format: sortable
    lib
    )

baci_add_dependency(lib_immersed ${LIBS_LIB_IMMERSED})

set(OBJS_ALG_IMMERSED_PROBLEMS
    # cmake-format: sortable
    immersed_problem_dyn.cpp
    immersed_problem_fsi_partitioned_immersed.cpp
    immersed_problem_immersed_base.cpp
    immersed_problem_immersed_partitioned_fsi_dirichletneumann.cpp
    )

baci_add_library(alg_immersed_problem ${OBJS_ALG_IMMERSED_PROBLEMS})

set(LIBS_ALG_IMMERSED_PROBLEMS
    # cmake-format: sortable
    alg_fsi
    )

baci_add_dependency(alg_immersed_problem ${LIBS_ALG_IMMERSED_PROBLEMS})
