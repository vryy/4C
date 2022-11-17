set(OBJS_DRT_TOPOPT
    topopt_algorithm.cpp
    topopt_dyn.cpp
    topopt_utils.cpp
    topopt_fluidAdjoint3_boundary.cpp
    topopt_fluidAdjoint3_impl.cpp
    topopt_fluidAdjoint3_impl_parameter.cpp
    topopt_optimizer_ele.cpp
    topopt_optimizer_ele_evaluate.cpp
    topopt_optimizer_ele_impl.cpp
    topopt_optimizer_ele_parameter.cpp
    topopt_fluidAdjointImplTimeIntegration.cpp
    topopt_fluidAdjointResulttest.cpp
    topopt_fluidAdjoint_timeint.cpp
    topopt_optimizer.cpp
    opti_GCMMA.cpp
    opti_resulttest.cpp
    )

baci_add_library(drt_opti ${OBJS_DRT_TOPOPT})

target_link_libraries(drt_opti PUBLIC drt_adapter)
