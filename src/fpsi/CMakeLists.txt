four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    ale
    config
    constraint
    core
    fluid
    fluid_ele
    fsi
    global_data
    inpar
    lib
    poroelast
    poroelast_scatra
    structure
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
