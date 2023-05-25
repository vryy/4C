baci_auto_define_module()

set(_fluid_turbulence_dependencies
    # cmake-format: sortable
    adapter
    core
    fluid
    fluid_ele
    inpar
    io
    lib
    linalg
    mat
    nurbs_discret
    scatra
    scatra_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_fluid_turbulence_dependencies})
