baci_auto_define_module()

set(_scatra_dependencies
    # cmake-format: sortable
    adapter
    binstrategy
    core
    fluid
    fluid_turbulence
    inpar
    io
    lib
    mat
    mortar
    nurbs_discret
    poromultiphase_scatra
    rebalance
    scatra_ele
    solver_nonlin_nox
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_scatra_dependencies})
