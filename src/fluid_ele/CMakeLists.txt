four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    bele
    config
    core
    cut
    fluid
    global_data
    immersed_problem
    inpar
    lib
    mat
    nurbs_discret
    poroelast
    xfem
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

# duplicate macros
four_c_unity_build_compile_separately(4C_fluid_ele_calc_poro.cpp)
four_c_unity_build_compile_separately(4C_fluid_ele_calc.cpp)
