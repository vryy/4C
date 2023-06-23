baci_auto_define_module()

set(_ssi_dependencies
    # cmake-format: sortable
    adapter
    contact
    core
    inpar
    io
    lib
    linalg
    mat
    rebalance
    scatra
    scatra_ele
    solver_nonlin_nox
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_ssi_dependencies})
