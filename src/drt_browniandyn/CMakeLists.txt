set(OBJS_DRT_LIB_BROWNIANDYN
    # cmake-format: sortable
    str_model_evaluator_browniandyn.cpp
    str_model_evaluator_browniandyndata.cpp
    )

baci_add_library(drt_lib_browniandyn ${OBJS_DRT_LIB_BROWNIANDYN})

set(LIBS_DRT_LIB_BROWNIANDYN
    # cmake-format: sortable
    alg_structure
    )

target_link_libraries(drt_lib_browniandyn PUBLIC ${LIBS_DRT_LIB_BROWNIANDYN})
