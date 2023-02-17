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

baci_add_library(levelset ${SOURCES})

set(LIBS_LEVELSET
    # cmake-format: sortable
    linalg
    )

baci_add_dependency(levelset ${LIBS_LEVELSET})
