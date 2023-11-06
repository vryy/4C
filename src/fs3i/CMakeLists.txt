baci_auto_define_module()

set(_fs3i_dependencies
    # cmake-format: sortable
    adapter
    ale
    core
    fluid
    fluid_ele
    fpsi
    fsi
    inpar
    io
    lib
    mat
    poroelast
    poroelast_scatra
    scatra
    scatra_ele
    ssi
    structure
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fs3i_dependencies})
