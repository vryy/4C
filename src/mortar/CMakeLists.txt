four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    contact
    core
    global_data
    inpar
    mat
    poroelast
    poroelast_scatra
    so3
    solver_nonlin_nox
    structure_new
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
