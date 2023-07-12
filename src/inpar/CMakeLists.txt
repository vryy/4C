baci_auto_define_module()

set(_inpar_dependencies
    # cmake-format: sortable
    contact_constitutivelaw
    core
    io
    lib
    solver_nonlin_nox
    tutorial
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_inpar_dependencies})
