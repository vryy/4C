set(SOURCES
    levelset_algorithm.cpp
    levelset_algorithm_utils.cpp
    levelset_algorithm_reinit.cpp
    levelset_dyn.cpp
    levelset_intersection_utils.cpp
    levelset_timint_ost.cpp
    levelset_timint_stat.cpp
    )

baci_add_library(drt_levelset ${SOURCES})
target_link_libraries(drt_levelset PUBLIC linalg)
