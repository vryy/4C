set(OBJS_LIB_BROWNIANDYN
    # cmake-format: sortable
    browniandyn_str_model_evaluator.cpp
    browniandyn_str_model_evaluator_data.cpp
    )

baci_add_library(lib_browniandyn ${OBJS_LIB_BROWNIANDYN})

set(LIBS_LIB_BROWNIANDYN
    # cmake-format: sortable
    alg_structure
    )

target_link_libraries(lib_browniandyn PUBLIC ${LIBS_LIB_BROWNIANDYN})
