baci_auto_define_module()

set(_poromultiphase_dependencies
    # cmake-format: sortable
    adapter
    inpar
    io
    lib
    linalg
    mat
    poroelast
    porofluidmultiphase
    porofluidmultiphase_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_poromultiphase_dependencies})
