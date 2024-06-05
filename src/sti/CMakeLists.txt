four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    core
    global_data
    inpar
    lib
    mat
    scatra
    scatra_ele
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
