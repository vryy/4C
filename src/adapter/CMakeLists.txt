four_c_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    ale
    beam3
    beaminteraction
    config
    constraint
    contact
    core
    fluid
    fluid_ele
    fluid_turbulence
    fluid_xfluid
    fpsi
    fsi
    global_data
    inpar
    levelset
    lib
    lubrication
    mat
    membrane
    mortar
    nurbs_discret
    pasi
    poroelast
    porofluidmultiphase
    scatra
    shell7p
    so3
    solid_3D_ele
    solver_nonlin_nox
    ssi
    structure
    structure_new
    thermo
    truss3
    w1
    xfem
    )

four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
