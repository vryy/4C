baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    bele
    config
    core
    cut
    fluid
    global_data
    immersed_problem
    inpar
    io
    lib
    mat
    nurbs_discret
    poroelast
    xfem
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
