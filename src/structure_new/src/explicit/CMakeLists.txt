baci_auto_define_module()

set(_structure_new_expl_dependencies
    # cmake-format: sortable
    solver_nonlin_nox
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_structure_new_expl_dependencies})
