set(OBJS_DRT_ELCH elch_algorithm.cpp elch_dyn.cpp elch_moving_boundary_algorithm.cpp)

baci_add_library(alg_elch ${OBJS_DRT_ELCH})

target_link_libraries(alg_elch PUBLIC linalg)
