baci_auto_define_module()

set(_geometry_pair_dependencies
    # cmake-format: sortable
    beam3
    core
    inpar
    lib
    nurbs_discret
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_geometry_pair_dependencies})
