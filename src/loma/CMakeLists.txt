four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    core
    fluid
    fluid_ele
    global_data
    inpar
    io
    lib
    scatra
    scatra_ele
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
