four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    contact
    core
    fluid
    fluid_ele
    fsi
    global_data
    inpar
    io
    lib
    mat
    mortar
    so3
    solid_poro_3D_ele
    structure
    w1
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
