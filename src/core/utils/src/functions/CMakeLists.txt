baci_auto_define_module()
set(LIBS_UTILS_FUNCTIONS
    # cmake-format: sortable
    core_utils_exceptions
    global_data
    io
    lib
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_UTILS_FUNCTIONS})
