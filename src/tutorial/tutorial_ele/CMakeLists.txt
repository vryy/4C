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

baci_add_dependency(alg_tutorial_ele ${LIBS_ALG_TUTORIAL_ELE})
