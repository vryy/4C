set(OBJS_SSTI
    # cmake-format: sortable
    ssti_algorithm.cpp
    ssti_dyn.cpp
    ssti_monolithic.cpp
    ssti_monolithic_assemble_strategy.cpp
    ssti_monolithic_evaluate_OffDiag.cpp
    ssti_resulttest.cpp
    ssti_utils.cpp
    )

baci_add_library(alg_ssti ${OBJS_SSTI})

set(LIBS_ALG_SSTI
    # cmake-format: sortable
    adapter
    alg_sti
    )

target_link_libraries(alg_ssti PUBLIC ${LIBS_ALG_SSTI})
