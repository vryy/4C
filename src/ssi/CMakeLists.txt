baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    contact
    core
    global_data
    inpar
    io
    lib
    mat
    scatra
    scatra_ele
    solver_nonlin_nox
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
