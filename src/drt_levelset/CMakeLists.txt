set(SOURCES
    # cmake-format: sortable
    levelset_algorithm.cpp
    levelset_algorithm_reinit.cpp
    levelset_algorithm_utils.cpp
    levelset_dyn.cpp
    levelset_intersection_utils.cpp
    levelset_timint_ost.cpp
    levelset_timint_stat.cpp
    )

baci_add_library(drt_levelset ${SOURCES})

set(LIBS_DRT_LEVELSET
    # cmake-format: sortable
    linalg
    )

target_link_libraries(drt_levelset PUBLIC ${LIBS_DRT_LEVELSET})
