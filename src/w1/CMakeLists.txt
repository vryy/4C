baci_auto_define_module()

set(_w1_dependencies
    # cmake-format: sortable
    contact
    core
    inpar
    lib
    linalg
    mat
    nurbs_discret
    poroelast
    so3
    structure_new
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_w1_dependencies})
