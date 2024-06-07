four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    constraint
    contact
    core
    fluid_xfluid
    fsi
    global_data
    inpar
    poroelast
    structure
    xfem
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
