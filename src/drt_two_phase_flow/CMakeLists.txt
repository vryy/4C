set(OBJS_DRT_TWO_PHASE two_phase_algorithm.cpp two_phase_dyn.cpp)
baci_add_library(alg_two_phase ${OBJS_DRT_TWO_PHASE})
target_link_libraries(alg_two_phase PUBLIC drt_lib alg_xfluid)
