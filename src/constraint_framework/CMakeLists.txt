baci_auto_define_module()

set(_constraint_framework_dependencies
    # cmake-format: sortable
    adapter
    core
    inpar
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_constraint_framework_dependencies})
