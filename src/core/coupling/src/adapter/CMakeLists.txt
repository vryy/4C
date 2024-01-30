baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    inpar
    io
    lib
    mortar
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
