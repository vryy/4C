baci_auto_define_module()

set(_ehl_dependencies
    # cmake-format: sortable
    adapter
    contact
    core
    inpar
    io
    lib
    linalg
    lubrication
    lubrication_ele
    mat
    mortar
    solver
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_ehl_dependencies})
