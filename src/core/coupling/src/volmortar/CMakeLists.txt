baci_auto_define_module()

set(LIBS_VOLMORTAR
    # cmake-format: sortable
    lib
    mortar
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_VOLMORTAR})
