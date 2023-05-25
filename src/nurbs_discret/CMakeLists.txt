baci_auto_define_module()

set(_nurbs_discret_dependencies
    # cmake-format: sortable
    core
    io
    lib
    linalg
    solver
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_nurbs_discret_dependencies})
