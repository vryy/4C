set(OBJS_STRUCT_ALE
    # cmake-format: sortable
    struct_ale_str_model_evaluator.cpp
    )

baci_add_library(lib_struct_ale ${OBJS_STRUCT_ALE})

set(LIBS_LIB_STRUCT_ALE
    # cmake-format: sortable
    alg_structure_new
    )

target_link_libraries(lib_struct_ale PUBLIC ${LIBS_LIB_STRUCT_ALE})
