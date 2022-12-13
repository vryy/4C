set(OBJS_DRT_ELEMAG
    # cmake-format: sortable
    elemag_dyn.cpp
    elemag_resulttest.cpp
    elemag_timeint.cpp
    elemag_utils_clonestrategy.cpp
    )

baci_add_library(alg_elemag ${OBJS_DRT_ELEMAG})

set(LIBS_ALG_ELEMEAG
    # cmake-format: sortable
    linalg
    )

target_link_libraries(alg_elemag PUBLIC ${LIBS_ALG_ELEMEAG})

set(OBJS_DRT_ELEMAG_LIB
    # cmake-format: sortable
    elemag_diff_ele.cpp
    elemag_diff_ele_calc.cpp
    elemag_ele.cpp
    elemag_ele_boundary_calc.cpp
    elemag_ele_calc.cpp
    elemag_ele_evaluate.cpp
    elemag_ele_factory.cpp
    elemag_ele_intfaces_calc.cpp
    )

baci_add_library(ele_elemag ${OBJS_DRT_ELEMAG_LIB})

set(LIBS_ELE_ELEMAG
    # cmake-format: sortable
    linalg
    )

target_link_libraries(ele_elemag PUBLIC ${LIBS_ELE_ELEMAG})
