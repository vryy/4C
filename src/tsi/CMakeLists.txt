baci_auto_define_module()

set(_tsi_dependencies
    # cmake-format: sortable
    adapter
    contact
    core
    inpar
    io
    lib
    linalg
    mortar
    so3
    solver
    structure_new
    thermo
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_tsi_dependencies})
