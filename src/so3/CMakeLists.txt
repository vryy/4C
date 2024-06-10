four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    contact
    core
    fluid_ele
    global_data
    immersed_problem
    inpar
    mat
    structure_new
    thermo
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

four_c_unity_build_compile_separately(4C_so3_scatra_evaluate.cpp)
