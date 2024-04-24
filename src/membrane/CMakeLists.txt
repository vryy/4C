baci_auto_define_module()

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

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

# duplicate template instantiations
baci_unity_build_compile_separately(4C_membrane_evaluate.cpp)
baci_unity_build_compile_separately(4C_membrane_input.cpp)
baci_unity_build_compile_separately(4C_membrane_line_evaluate.cpp)
baci_unity_build_compile_separately(4C_membrane_line_input.cpp)
baci_unity_build_compile_separately(4C_membrane_scatra_evaluate.cpp)
