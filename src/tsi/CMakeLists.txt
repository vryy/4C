baci_auto_define_module()

set(_tsi_dependencies
    # cmake-format: sortable
    adapter
    contact
    core
    inpar
    io
    lib
    mortar
    so3
    structure_new
    thermo
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_tsi_dependencies})
