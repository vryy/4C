four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    global_data
    inpar
    io
    mat
    matelast
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

four_c_unity_build_compile_separately(4C_mixture_constituent_remodelfiber_expl.cpp)
four_c_unity_build_compile_separately(4C_mixture_constituent_remodelfiber_impl.cpp)
four_c_unity_build_compile_separately(4C_mixture_constituent_full_constrained_mixture_fiber.cpp)
