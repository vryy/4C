set(OBJS_VOLMORTAR
    # cmake-format: sortable
    coupling_volmortar.cpp
    coupling_volmortar_cell.cpp
    coupling_volmortar_integrator.cpp
    coupling_volmortar_utils.cpp
    )

set(LIBS_VOLMORTAR
    # cmake-format: sortable
    lib
    lib_mortar
    mortar
    )

baci_add_library_to_module(
  MODULE
  coupling
  NAME
  volmortar
  SOURCES
  ${OBJS_VOLMORTAR}
  DEPENDENCIES
  ${LIBS_VOLMORTAR}
  )
