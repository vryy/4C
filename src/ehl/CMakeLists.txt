baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    contact
    core
    global_data
    inpar
    io
    lib
    lubrication
    lubrication_ele
    mat
    mortar
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
