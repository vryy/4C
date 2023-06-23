baci_auto_define_module()

set(_immersed_problem_dependencies
    # cmake-format: sortable
    adapter
    ale
    core
    cut
    fluid_ele
    fsi
    inpar
    io
    lib
    linalg
    mortar
    structure
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_immersed_problem_dependencies})
