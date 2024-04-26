four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    global_data
    inpar
    io
    lib
    mat
    so3
    structure_new
    thermo
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

# duplicate template instantiations
four_c_unity_build_compile_separately(4C_membrane_evaluate.cpp)
four_c_unity_build_compile_separately(4C_membrane_input.cpp)
four_c_unity_build_compile_separately(4C_membrane_line_evaluate.cpp)
four_c_unity_build_compile_separately(4C_membrane_line_input.cpp)
four_c_unity_build_compile_separately(4C_membrane_scatra_evaluate.cpp)
