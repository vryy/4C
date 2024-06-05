four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    core
    fluid_ele
    global_data
    inpar
    lib
    mat
    poroelast
    scatra
    scatra_ele
    so3
    solid_poro_3D_ele
    w1
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
