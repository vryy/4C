four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    ale
    config
    core
    fluid
    fluid_ele
    fpsi
    fsi
    global_data
    inpar
    mat
    poroelast
    poroelast_scatra
    scatra
    scatra_ele
    ssi
    structure
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
