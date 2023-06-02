set(OBJS_PRE_LOCSYS
    # cmake-format: sortable
    pre_locsys.cpp
    )

add_executable(pre_locsys EXCLUDE_FROM_ALL ${OBJS_PRE_LOCSYS})
set_target_properties(pre_locsys PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_target_requires_cyclic_legacy_linking(pre_locsys)
