baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    contact
    core
    fluid_ele
    global_data
    inpar
    io
    lib
    mat
    nurbs_discret
    poroelast
    so3
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

baci_unity_build_compile_separately(4C_w1_poro_evaluate.cpp)
baci_unity_build_compile_separately(4C_w1_poro_p1_evaluate.cpp)
