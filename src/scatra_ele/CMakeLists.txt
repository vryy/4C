baci_auto_define_module()

set(_scatra_ele_dependencies
    # cmake-format: sortable
    core
    fiber
    fluid
    fluid_ele
    inpar
    lib
    mat
    nurbs_discret
    porofluidmultiphase_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_scatra_ele_dependencies})
