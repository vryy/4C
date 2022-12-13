set(OBJS_DRT_STI
    # cmake-format: sortable
    sti_algorithm.cpp
    sti_clonestrategy.cpp
    sti_dyn.cpp
    sti_monolithic.cpp
    sti_monolithic_evaluate_OffDiag.cpp
    sti_partitioned.cpp
    sti_resulttest.cpp
    )

baci_add_library(alg_sti ${OBJS_DRT_STI})

set(LIBS_ALG_STI
    # cmake-format: sortable
    drt_adapter
    drt_lib
    )

target_link_libraries(alg_sti PUBLIC ${LIBS_ALG_STI})
