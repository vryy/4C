baci_auto_define_module()

set(_structure_new_output_dependencies
    # cmake-format: sortable
    io
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_structure_new_output_dependencies})
