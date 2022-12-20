set(OBJS_TUTORIAL_ELE
    # cmake-format: sortable
    tutorial_ele.cpp
    tutorial_ele_calc.cpp
    )

baci_add_library(alg_tutorial_ele ${OBJS_TUTORIAL_ELE})

set(LIBS_ALG_TUTORIAL_ELE
    # cmake-format: sortable
    alg_tutorial
    )

target_link_libraries(alg_tutorial_ele PUBLIC ${LIBS_ALG_TUTORIAL_ELE})
