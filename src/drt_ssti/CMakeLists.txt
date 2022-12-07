set(OBJS_DRT_SSTI
    ssti_algorithm.cpp
    ssti_dyn.cpp
    ssti_monolithic.cpp
    ssti_monolithic_assemble_strategy.cpp
    ssti_monolithic_evaluate_OffDiag.cpp
    ssti_resulttest.cpp
    ssti_utils.cpp
    )

baci_add_library(alg_ssti ${OBJS_DRT_SSTI})
target_link_libraries(alg_ssti PUBLIC drt_adapter)
