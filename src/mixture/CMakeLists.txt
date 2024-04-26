baci_auto_define_module()

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

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

baci_unity_build_compile_separately(4C_mixture_constituent_remodelfiber_expl.cpp)
baci_unity_build_compile_separately(4C_mixture_constituent_remodelfiber_impl.cpp)
baci_unity_build_compile_separately(4C_mixture_constituent_full_constrained_mixture_fiber.cpp)
