four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    contact
    core
    fluid_ele
    global_data
    inpar
    mat
    poroelast
    so3
    structure_new
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

four_c_unity_build_compile_separately(4C_w1_poro_evaluate.cpp)
four_c_unity_build_compile_separately(4C_w1_poro_p1_evaluate.cpp)
