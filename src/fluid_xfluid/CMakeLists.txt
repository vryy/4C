baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    core
    cut
    fluid
    fluid_ele
    global_data
    inpar
    io
    levelset
    lib
    mat
    xfem
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
