baci_auto_define_module()

set(_io_dependencies
    # cmake-format: sortable
    beam3
    config
    core
    fluid_ele
    geometry
    inpar
    lib
    linalg
    nurbs_discret
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_io_dependencies})
