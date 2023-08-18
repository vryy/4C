baci_auto_define_module()

set(_ale_dependencies
    # cmake-format: sortable
    adapter
    contact
    core
    fluid_ele
    inpar
    io
    lib
    mat
    mortar
    nurbs_discret
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_ale_dependencies})
