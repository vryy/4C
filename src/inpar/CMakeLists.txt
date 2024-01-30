baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    contact_constitutivelaw
    core
    global_data
    io
    lib
    solver_nonlin_nox
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
