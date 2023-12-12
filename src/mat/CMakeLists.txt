baci_auto_define_module()

set(_mat_dependencies
    # cmake-format: sortable
    core
    fiber
    inpar
    io
    lib
    matelast
    mixture
    porofluidmultiphase_ele
    red_airways
    scatra
    scatra_ele
    so3
    stru_multi
    structure_new
    tsi
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_mat_dependencies})
