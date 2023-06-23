baci_auto_define_module()

set(_poromultiphase_scatra_dependencies
    # cmake-format: sortable
    adapter
    art_net
    core
    inpar
    io
    lib
    linalg
    mat
    poroelast
    poroelast_scatra
    porofluidmultiphase
    porofluidmultiphase_ele
    poromultiphase
    scatra
    scatra_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_poromultiphase_scatra_dependencies})
