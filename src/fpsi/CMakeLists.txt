baci_auto_define_module()

set(_fpsi_dependencies
    # cmake-format: sortable
    adapter
    ale
    constraint
    core
    fluid
    fluid_ele
    fsi
    inpar
    io
    lib
    linalg
    poroelast
    poroelast_scatra
    structure
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fpsi_dependencies})
