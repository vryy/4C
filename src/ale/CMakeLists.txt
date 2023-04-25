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
    adapter
    ele_ale
    inpar
    )

baci_add_dependency(alg_ale ${LIBS_ALG_ALE})

set(OBJS_ELE_ALE
    # cmake-format: sortable
    ale_ale2.cpp
    ale_ale2_evaluate.cpp
    ale_ale2_input.cpp
    ale_ale2_line.cpp
    ale_ale2_line_evaluate.cpp
    ale_ale2_nurbs.cpp
    ale_ale3.cpp
    ale_ale3_evaluate.cpp
    ale_ale3_input.cpp
    ale_ale3_nurbs.cpp
    ale_ale3_surface.cpp
    ale_ale3_surface_evaluate.cpp
    )

baci_add_library(ele_ale ${OBJS_ELE_ALE})

set(LIBS_ELE_ALE
    # cmake-format: sortable
    ele_so
    lib
    mat
    nurbs_discret
    )

baci_add_dependency(ele_ale ${LIBS_ELE_ALE})
