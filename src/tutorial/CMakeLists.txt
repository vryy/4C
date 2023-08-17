baci_auto_define_module()

set(_tutorial_dependencies
    # cmake-format: sortable
    core
    inpar
    lib
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_tutorial_dependencies})
