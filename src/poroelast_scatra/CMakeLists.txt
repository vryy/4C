baci_auto_define_module()

set(_poroelast_scatra_dependencies
    # cmake-format: sortable
    adapter
    core
    fluid_ele
    inpar
    io
    lib
    mat
    poroelast
    scatra
    scatra_ele
    so3
    solid_poro_ele
    w1
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_poroelast_scatra_dependencies})
