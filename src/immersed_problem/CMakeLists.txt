baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    ale
    config
    core
    cut
    fluid_ele
    fsi
    global_data
    inpar
    io
    lib
    mortar
    structure
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
