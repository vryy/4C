set(OBJS_DRT_LUBRICATION
    lubrication_dyn.cpp
    lubrication_resulttest.cpp
    lubrication_timint_implicit.cpp
    lubrication_timint_stat.cpp
    )

baci_add_library(alg_lubrication ${OBJS_DRT_LUBRICATION})
target_link_libraries(alg_lubrication PUBLIC drt_lib_inpar drt_adapter ele_lubrication)
