set(OBJS_DRT_LOMA loma_algorithm.cpp loma_dyn.cpp)

baci_add_library(alg_loma ${OBJS_DRT_LOMA})
target_link_libraries(alg_loma PUBLIC linalg)
