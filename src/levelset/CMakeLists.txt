baci_auto_define_module()

set(_levelset_dependencies
    # cmake-format: sortable
    adapter
    core
    cut
    inpar
    io
    lib
    linalg
    scatra
    scatra_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_levelset_dependencies})
