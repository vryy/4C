baci_auto_define_module()

set(_ssti_dependencies
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
    so3
    ssi
    sti
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_ssti_dependencies})
