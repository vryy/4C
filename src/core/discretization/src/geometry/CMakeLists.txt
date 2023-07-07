baci_auto_define_module()

set(_geometry_dependencies
    # cmake-format: sortable
    cut
    io
    lib
    linalg
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_geometry_dependencies})
