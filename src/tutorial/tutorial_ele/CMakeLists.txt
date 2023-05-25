baci_auto_define_module()

set(LIBS
    # cmake-format: sortable
    tutorial
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS})
