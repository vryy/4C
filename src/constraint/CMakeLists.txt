baci_auto_define_module()

set(_constraint_dependencies
    # cmake-format: sortable
    adapter
    contact
    core
    inpar
    io
    lib
    linalg
    solver_nonlin_nox
    structure_new
    truss3
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_constraint_dependencies})
