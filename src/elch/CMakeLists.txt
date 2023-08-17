baci_auto_define_module()

set(_elch_dependencies
    # cmake-format: sortable
    adapter
    ale
    fluid
    inpar
    io
    lib
    scatra
    scatra_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_elch_dependencies})
