baci_auto_define_module()

set(_shell7p_dependencies
    # cmake-format: sortable
    core
    inpar
    io
    lib
    mat
    matelast
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_shell7p_dependencies})
