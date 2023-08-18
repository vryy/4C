baci_auto_define_module()

set(_geometry_dependencies
    # cmake-format: sortable
    core_linalg
    cut
    io
    lib
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_geometry_dependencies})
