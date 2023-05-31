baci_auto_define_module()

set(_sti_dependencies
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
    solver
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_sti_dependencies})
