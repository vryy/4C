baci_auto_define_module()

set(_art_net_dependencies
    # cmake-format: sortable
    adapter
    core
    inpar
    io
    lib
    linalg
    mat
    scatra
    scatra_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_art_net_dependencies})
