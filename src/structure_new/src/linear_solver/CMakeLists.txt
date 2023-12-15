baci_auto_define_module()

set(_structure_new_lin_solver_dependencies
    # cmake-format: sortable
    beam3
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_structure_new_lin_solver_dependencies})
