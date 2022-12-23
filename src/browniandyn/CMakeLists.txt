set(OBJS_LIB_BROWNIANDYN
    # cmake-format: sortable
    str_model_evaluator_browniandyn.cpp
    str_model_evaluator_browniandyndata.cpp
    )

baci_add_library(lib_browniandyn ${OBJS_LIB_BROWNIANDYN})

set(LIBS_LIB_BROWNIANDYN
    # cmake-format: sortable
    alg_structure
    )

target_link_libraries(lib_browniandyn PUBLIC ${LIBS_LIB_BROWNIANDYN})
