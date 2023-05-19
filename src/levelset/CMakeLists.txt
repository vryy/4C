baci_auto_define_module()

set(_levelset_dependencies
    # cmake-format: sortable
    adapter
    cut
    geometry
    inpar
    io
    lib
    linalg
    scatra
    scatra_ele
    solver
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_levelset_dependencies})
