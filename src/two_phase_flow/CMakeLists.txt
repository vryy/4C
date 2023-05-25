baci_auto_define_module()

set(_two_phase_flow_dependencies
    # cmake-format: sortable
    adapter
    fluid
    fluid_xfluid
    inpar
    levelset
    lib
    linalg
    scatra
    scatra_ele
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_two_phase_flow_dependencies})
