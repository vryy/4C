baci_auto_define_module()

set(_poroelast_dependencies
    # cmake-format: sortable
    adapter
    contact
    core
    fluid
    fluid_ele
    fsi
    inpar
    io
    lib
    mat
    mortar
    so3
    structure
    w1
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_poroelast_dependencies})
