set(OBJS_DRT_STI
    sti_algorithm.cpp
    sti_clonestrategy.cpp
    sti_dyn.cpp
    sti_monolithic.cpp
    sti_monolithic_evaluate_OffDiag.cpp
    sti_partitioned.cpp
    sti_resulttest.cpp
    )

baci_add_library(alg_sti ${OBJS_DRT_STI})
target_link_libraries(alg_sti PUBLIC drt_adapter drt_lib)
