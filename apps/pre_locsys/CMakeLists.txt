set(OBJS_PRE_LOCSYS
    # cmake-format: sortable
    baci_pre_locsys.cpp
    )

add_executable(pre_locsys EXCLUDE_FROM_ALL ${OBJS_PRE_LOCSYS})
set_target_properties(pre_locsys PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_set_up_executable(pre_locsys)
