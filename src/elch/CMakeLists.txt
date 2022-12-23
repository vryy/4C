set(OBJS_ELCH
    # cmake-format: sortable
    elch_algorithm.cpp
    elch_dyn.cpp
    elch_moving_boundary_algorithm.cpp
    )

baci_add_library(alg_elch ${OBJS_ELCH})

set(LIBS_ALG_ELCH
    # cmake-format: sortable
    linalg
    )

target_link_libraries(alg_elch PUBLIC ${LIBS_ALG_ELCH})
