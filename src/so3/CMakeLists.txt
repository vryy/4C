baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    contact
    core
    fiber
    fluid_ele
    global_data
    immersed_problem
    inpar
    io
    lib
    mat
    nurbs_discret
    structure_new
    thermo
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
