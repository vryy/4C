baci_auto_define_module()

set(LIBS_ADAPTER_COUPLING
    # cmake-format: sortable
    comm
    constraint
    core_coupling_volmortar
    io
    levelset
    lib
    linalg
    mortar
    nurbs_discret
    )

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${LIBS_ADAPTER_COUPLING})
