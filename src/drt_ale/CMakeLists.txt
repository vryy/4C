set(OBJS_ALG_ALE
    # cmake-format: sortable
    ale.cpp
    ale_dyn.cpp
    ale_meshsliding.cpp
    ale_meshtying.cpp
    ale_resulttest.cpp
    ale_utils_clonestrategy.cpp
    ale_utils_mapextractor.cpp
    )

baci_add_library(alg_ale ${OBJS_ALG_ALE})

set(LIBS_ALG_ALE
    # cmake-format: sortable
    drt_adapter
    drt_inpar
    ele_ale
    )

target_link_libraries(alg_ale PUBLIC ${LIBS_ALG_ALE})

set(OBJS_ELE_ALE
    # cmake-format: sortable
    ale2.cpp
    ale2_evaluate.cpp
    ale2_input.cpp
    ale2_line.cpp
    ale2_line_evaluate.cpp
    ale2_nurbs.cpp
    ale3.cpp
    ale3_evaluate.cpp
    ale3_input.cpp
    ale3_nurbs.cpp
    ale3_surface.cpp
    ale3_surface_evaluate.cpp
    )

baci_add_library(ele_ale ${OBJS_ELE_ALE})

set(LIBS_ELE_ALE
    # cmake-format: sortable
    drt_lib
    drt_mat
    drt_nurbs_discret
    ele_so
    )

target_link_libraries(ele_ale PUBLIC ${LIBS_ELE_ALE})
