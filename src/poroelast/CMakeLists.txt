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
    linalg
    mat
    mortar
    so3
    solver
    structure
    w1
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_poroelast_dependencies})
