set(OBJS_TOPOPT
    # cmake-format: sortable
    opti_GCMMA.cpp
    opti_resulttest.cpp
    opti_topopt_algorithm.cpp
    opti_topopt_dyn.cpp
    opti_topopt_fluidAdjoint3_boundary.cpp
    opti_topopt_fluidAdjoint3_impl.cpp
    opti_topopt_fluidAdjoint3_impl_parameter.cpp
    opti_topopt_fluidAdjoint_timeint.cpp
    opti_topopt_fluidAdjointImplTimeIntegration.cpp
    opti_topopt_fluidAdjointResulttest.cpp
    opti_topopt_optimizer.cpp
    opti_topopt_optimizer_ele.cpp
    opti_topopt_optimizer_ele_evaluate.cpp
    opti_topopt_optimizer_ele_impl.cpp
    opti_topopt_optimizer_ele_parameter.cpp
    opti_topopt_utils.cpp
    )

baci_add_library(opti ${OBJS_TOPOPT})

set(LIBS_OPTI
    # cmake-format: sortable
    adapter
    fem_general
    inpar
    lib
    )

target_link_libraries(opti PUBLIC ${LIBS_OPTI})
