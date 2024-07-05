four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    cardiovascular0d
    config
    constraint
    contact
    core
    inpar
    structure_new
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

four_c_unity_build_compile_separately(4C_solver_nonlin_nox_inner_statustest_filter.cpp)
