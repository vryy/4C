baci_auto_define_module()

set(LIBS_SEARCHING
    # cmake-format: sortable
    lib
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_SEARCHING})
