set(OBJS_DRT_STRUCTURE_XSTRUCTURE
    xstr_xstructure_state.cpp
    xstr_xstructure_structure_state.cpp
    xstr_multi_discretization_wrapper.cpp
    xstr_multi_io.cpp
    str_timint_base_xfem.cpp
    str_timint_implicitbase_xfem.cpp
    )

baci_add_library(alg_xstructure ${OBJS_DRT_STRUCTURE_XSTRUCTURE})
target_link_libraries(alg_xstructure PUBLIC alg_structure)
