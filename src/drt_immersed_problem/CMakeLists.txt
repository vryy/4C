set(OBJS_DRT_IMMERSED_LIB immersed_node.cpp)

baci_add_library(drt_lib_immersed ${OBJS_DRT_IMMERSED_LIB})
target_link_libraries(drt_lib_immersed PUBLIC drt_lib)

set(OBJS_DRT_IMMERSED_PROBLEM
    immersed_base.cpp
    immersed_problem_dyn.cpp
    immersed_partitioned_fsi_dirichletneumann.cpp
    fsi_partitioned_immersed.cpp
    )

baci_add_library(alg_immersed_problem ${OBJS_DRT_IMMERSED_PROBLEM})
target_link_libraries(alg_immersed_problem PUBLIC alg_fsi)
