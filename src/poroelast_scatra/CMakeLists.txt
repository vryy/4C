baci_auto_define_module()

set(_poroelast_scatra_dependencies
    # cmake-format: sortable
    adapter
    core
    fluid_ele
    inpar
    io
    lib
    linalg
    mat
    poroelast
    scatra
    scatra_ele
    so3
    w1
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_poroelast_scatra_dependencies})
