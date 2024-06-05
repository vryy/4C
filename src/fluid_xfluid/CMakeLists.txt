four_c_auto_define_module()

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
    levelset
    lib
    mat
    xfem
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
