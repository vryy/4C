baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    fiber
    global_data
    inpar
    io
    lib
    mat
    nurbs_discret
    so3
    structure_new
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

baci_unity_build_compile_separately(4C_solid_3D_ele_calc_eas.cpp)
