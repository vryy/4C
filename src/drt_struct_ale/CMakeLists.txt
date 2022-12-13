set(OBJS_DRT_STRUCT_ALE
    # cmake-format: sortable
    struct_ale_str_model_evaluator.cpp
    )

baci_add_library(drt_lib_struct_ale ${OBJS_DRT_STRUCT_ALE})

set(LIBS_DRT_LIB_STRUCT_ALE
    # cmake-format: sortable
    alg_structure_new
    )

target_link_libraries(drt_lib_struct_ale PUBLIC ${LIBS_DRT_LIB_STRUCT_ALE})
