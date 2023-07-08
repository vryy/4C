baci_auto_define_module()

set(_mixture_dependencies
    # cmake-format: sortable
    inpar
    lib
    mat
    matelast
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_mixture_dependencies})
