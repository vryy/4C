four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    global_data
    inpar
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

# Clashing definitions of `View` inside Trilinos
four_c_unity_build_compile_separately(4C_linear_solver_method_direct.cpp)
