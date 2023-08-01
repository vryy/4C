set(OBJS_POST_S8CONVERT
    # cmake-format: sortable
    baci_post_s8convert.cpp
    )

add_executable(post_s8convert EXCLUDE_FROM_ALL ${OBJS_POST_S8CONVERT})
set_target_properties(post_s8convert PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_add_dependency(post_s8convert post_common)
baci_target_requires_cyclic_legacy_linking(post_s8convert)
