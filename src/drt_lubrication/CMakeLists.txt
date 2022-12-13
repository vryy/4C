set(OBJS_DRT_LUBRICATION
    # cmake-format: sortable
    lubrication_dyn.cpp
    lubrication_resulttest.cpp
    lubrication_timint_implicit.cpp
    lubrication_timint_stat.cpp
    )

baci_add_library(alg_lubrication ${OBJS_DRT_LUBRICATION})

set(LIBS_ALG_LUBRICATION
    # cmake-format: sortable
    drt_adapter
    drt_inpar
    ele_lubrication
    )

target_link_libraries(alg_lubrication PUBLIC ${LIBS_ALG_LUBRICATION})
