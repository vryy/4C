baci_auto_define_module()

set(_LIBS
    # cmake-format: sortable
    core_linalg
    core_utils_exceptions
    io
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_LIBS})
