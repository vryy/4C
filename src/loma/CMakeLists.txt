baci_auto_define_module()

set(_loma_dependencies
    # cmake-format: sortable
    adapter
    fluid
    fluid_ele
    inpar
    io
    lib
    linalg
    scatra
    scatra_ele
    solver
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_loma_dependencies})
