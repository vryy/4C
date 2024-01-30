baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    cardiovascular0d
    config
    constraint
    contact
    contact_aug
    core
    inpar
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
