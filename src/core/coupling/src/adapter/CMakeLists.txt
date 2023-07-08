baci_auto_define_module()

set(LIBS_ADAPTER_COUPLING
    # cmake-format: sortable
    comm
    constraint
    core_coupling_volmortar
    core_linalg
    io
    levelset
    lib
    mortar
    nurbs_discret
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_ADAPTER_COUPLING})
