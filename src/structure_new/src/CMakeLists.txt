baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    beam3
    beaminteraction
    config
    contact
    core
    global_data
    inpar
    io
    lib
    solver_nonlin_nox
    structure_new
    timestepping
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

baci_unity_build_compile_separately(4C_structure_new_timint_basedataglobalstate.cpp)
