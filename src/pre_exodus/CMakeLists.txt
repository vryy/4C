set(OBJS_PRE_EXODUS
    # cmake-format: sortable
    baci_pre_exodus.cpp
    baci_pre_exodus_centerline.cpp
    baci_pre_exodus_readbc.cpp
    baci_pre_exodus_reader.cpp
    baci_pre_exodus_soshextrusion.cpp
    baci_pre_exodus_validate.cpp
    baci_pre_exodus_writedat.cpp
    )

add_executable(pre_exodus EXCLUDE_FROM_ALL ${OBJS_PRE_EXODUS})
set_target_properties(pre_exodus PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} post_common)
baci_target_requires_cyclic_legacy_linking(pre_exodus)
