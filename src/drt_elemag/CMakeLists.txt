set(OBJS_DRT_ELEMAG
    elemag_dyn.cpp elemag_timeint.cpp elemag_resulttest.cpp elemag_utils_clonestrategy.cpp
    )

baci_add_library(alg_elemag ${OBJS_DRT_ELEMAG})

target_link_libraries(alg_elemag PUBLIC linalg)

set(OBJS_DRT_ELEMAG_LIB
    elemag_ele_boundary_calc.cpp
    elemag_ele_calc.cpp
    elemag_diff_ele_calc.cpp
    elemag_ele_evaluate.cpp
    elemag_ele_factory.cpp
    elemag_ele_intfaces_calc.cpp
    elemag_ele.cpp
    elemag_diff_ele.cpp
    )

baci_add_library(ele_elemag ${OBJS_DRT_ELEMAG_LIB})

target_link_libraries(ele_elemag PUBLIC linalg)
