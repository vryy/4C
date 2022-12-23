set(OBJS_TOPOPT
    # cmake-format: sortable
    opti_GCMMA.cpp
    opti_resulttest.cpp
    topopt_algorithm.cpp
    topopt_dyn.cpp
    topopt_fluidAdjoint3_boundary.cpp
    topopt_fluidAdjoint3_impl.cpp
    topopt_fluidAdjoint3_impl_parameter.cpp
    topopt_fluidAdjoint_timeint.cpp
    topopt_fluidAdjointImplTimeIntegration.cpp
    topopt_fluidAdjointResulttest.cpp
    topopt_optimizer.cpp
    topopt_optimizer_ele.cpp
    topopt_optimizer_ele_evaluate.cpp
    topopt_optimizer_ele_impl.cpp
    topopt_optimizer_ele_parameter.cpp
    topopt_utils.cpp
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
