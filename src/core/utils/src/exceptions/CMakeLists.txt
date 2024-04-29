four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
four_c_add_external_dependency(${AUTO_DEFINED_MODULE_NAME} ${CMAKE_DL_LIBS})
