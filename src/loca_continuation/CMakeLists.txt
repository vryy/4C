baci_auto_define_module()

set(_loca_continuation_dependencies
    # cmake-format: sortable
    lib
    solver_nonlin_nox
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_loca_continuation_dependencies})
