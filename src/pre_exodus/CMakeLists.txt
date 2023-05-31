set(OBJS_PRE_EXODUS
    # cmake-format: sortable
    pre_exodus.cpp
    pre_exodus_centerline.cpp
    pre_exodus_readbc.cpp
    pre_exodus_reader.cpp
    pre_exodus_soshextrusion.cpp
    pre_exodus_validate.cpp
    pre_exodus_writedat.cpp
    )

add_executable(pre_exodus EXCLUDE_FROM_ALL ${OBJS_PRE_EXODUS})
set_target_properties(pre_exodus PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} post_common)
baci_link_executable_with_all(pre_exodus)
