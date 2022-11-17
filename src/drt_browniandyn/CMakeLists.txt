set(OBJS_DRT_LIB_BROWNIANDYN
    str_model_evaluator_browniandyn.cpp str_model_evaluator_browniandyndata.cpp
    )

baci_add_library(drt_lib_browniandyn ${OBJS_DRT_LIB_BROWNIANDYN})
target_link_libraries(drt_lib_browniandyn PUBLIC alg_structure)
