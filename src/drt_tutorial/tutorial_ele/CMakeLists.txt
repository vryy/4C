set(OBJS_DRT_TUTORIAL_ELE tutorial_ele_calc.cpp tutorial_ele.cpp)

baci_add_library(alg_tutorial_ele ${OBJS_DRT_TUTORIAL_ELE})
target_link_libraries(alg_tutorial_ele PUBLIC alg_tutorial)
