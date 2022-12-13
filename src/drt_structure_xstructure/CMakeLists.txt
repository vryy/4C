set(OBJS_DRT_STRUCTURE_XSTRUCTURE
    # cmake-format: sortable
    str_timint_base_xfem.cpp
    str_timint_implicitbase_xfem.cpp
    xstr_multi_discretization_wrapper.cpp
    xstr_multi_io.cpp
    xstr_xstructure_state.cpp
    xstr_xstructure_structure_state.cpp
    )

baci_add_library(alg_xstructure ${OBJS_DRT_STRUCTURE_XSTRUCTURE})

set(LIBS_ALG_XSTRUCTURE
    # cmake-format: sortable
    alg_structure
    )

target_link_libraries(alg_xstructure PUBLIC ${LIBS_ALG_XSTRUCTURE})
