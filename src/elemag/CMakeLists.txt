baci_auto_define_module()

set(_elemag_dependencies
    # cmake-format: sortable
    adapter
    comm
    core
    cut
    inpar
    io
    lib
    linalg
    mat
    scatra
    scatra_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_elemag_dependencies})
