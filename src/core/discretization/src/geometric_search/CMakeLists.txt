baci_auto_define_module()

set(LIBS_GEOMETRIC_SEARCH
    # cmake-format: sortable
    io
    lib
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_GEOMETRIC_SEARCH})
