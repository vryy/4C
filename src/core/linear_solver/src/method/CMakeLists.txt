baci_auto_define_module()

set(LIBS_LINEAR_SOLVER_METHOD
    # cmake-format: sortable
    linalg
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_LINEAR_SOLVER_METHOD})
