baci_auto_define_module()

set(_membrane_dependencies
    # cmake-format: sortable
    core
    inpar
    lib
    mat
    structure_new
    thermo
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_membrane_dependencies})
