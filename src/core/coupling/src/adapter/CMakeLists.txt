set(OBJS_ADAPTER_COUPLING
    # cmake-format: sortable
    coupling_adapter.cpp
    coupling_adapter_converter.cpp
    coupling_adapter_mortar.cpp
    coupling_adapter_volmortar.cpp
    )

set(LIBS_ADAPTER_COUPLING
    # cmake-format: sortable
    comm
    constraint
    coupling_volmortar
    io
    levelset
    lib
    lib_mortar
    linalg
    nurbs_discret
    solver
    )

baci_add_library_to_module(
  MODULE
  coupling
  NAME
  adapter_coupling
  SOURCES
  ${OBJS_ADAPTER_COUPLING}
  DEPENDENCIES
  ${LIBS_ADAPTER_COUPLING}
  )
