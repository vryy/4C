set(OBJS_STRUCTURE_XSTRUCTURE
    # cmake-format: sortable
    structure_xstructure_state.cpp
    structure_xstructure_str_timint_base_xfem.cpp
    structure_xstructure_str_timint_implicitbase_xfem.cpp
    structure_xstructure_xstr_multi_discretization_wrapper.cpp
    structure_xstructure_xstr_multi_io.cpp
    structure_xstructure_xstr_state.cpp
    )

baci_add_library(alg_xstructure ${OBJS_STRUCTURE_XSTRUCTURE})

set(LIBS_ALG_XSTRUCTURE
    # cmake-format: sortable
    alg_structure
    )

target_link_libraries(alg_xstructure PUBLIC ${LIBS_ALG_XSTRUCTURE})
