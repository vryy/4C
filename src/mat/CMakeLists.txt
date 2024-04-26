four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    fiber
    global_data
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

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
