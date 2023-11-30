baci_auto_define_module()

set(_structure_new_model_eval_dependencies
    # cmake-format: sortable
    browniandyn
    constraint
    constraint_framework
    inpar
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_structure_new_model_eval_dependencies})
