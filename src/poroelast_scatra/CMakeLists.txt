baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    config
    core
    fluid_ele
    global_data
    inpar
    io
    lib
    mat
    poroelast
    scatra
    scatra_ele
    so3
    solid_poro_ele
    w1
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
