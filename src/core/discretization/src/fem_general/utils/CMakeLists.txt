baci_auto_define_module()

set(LIBS_FEM_GENERAL_UTILS
    # cmake-format: sortable
    cut
    geometry
    lib
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_FEM_GENERAL_UTILS})
