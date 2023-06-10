baci_auto_define_module()

set(_cardiovascular0d_dependencies
    # cmake-format: sortable
    adapter
    core
    inpar
    io
    lib
    linalg
    mor
    so3
    solver_nonlin_nox
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_cardiovascular0d_dependencies})
