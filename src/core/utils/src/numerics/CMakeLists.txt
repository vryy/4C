baci_auto_define_module()
set(LIBS_UTILS_NUMERICS
    # cmake-format: sortable
    core_linalg
    core_utils_exceptions
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_UTILS_NUMERICS})
