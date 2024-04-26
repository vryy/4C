baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    config
    core
    global_data
    inpar
    io
    lib
    mat
    nurbs_discret
    so3
    structure_new
    thermo
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})

baci_unity_build_compile_separately(4C_so3_plast_ssn_eas.cpp)
baci_unity_build_compile_separately(4C_so3_plast_ssn_evaluate.cpp)
