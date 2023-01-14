set(OBJS_TWO_PHASE
    # cmake-format: sortable
    two_phase_flow_algorithm.cpp
    two_phase_flow_dyn.cpp
    )

baci_add_library(alg_two_phase ${OBJS_TWO_PHASE})

set(LIBS_ALG_TWO_PHASE
    # cmake-format: sortable
    alg_xfluid
    lib
    )

target_link_libraries(alg_two_phase PUBLIC ${LIBS_ALG_TWO_PHASE})
