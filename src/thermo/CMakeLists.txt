set(OBJS_THERMO_ELE
    # cmake-format: sortable
    thermo_ele_boundary_impl.cpp
    thermo_ele_impl.cpp
    thermo_element.cpp
    thermo_element_boundary_evaluate.cpp
    thermo_element_evaluate.cpp
    thermo_element_input.cpp
    )

baci_add_library(ele_thermo ${OBJS_THERMO_ELE})

set(LIBS_ELE_THERMO
    # cmake-format: sortable
    inpar
    lib
    )

baci_add_dependency(ele_thermo ${LIBS_ELE_THERMO})

set(OBJS_THERMO
    # cmake-format: sortable
    thermo_aux.cpp
    thermo_dyn.cpp
    thermo_resulttest.cpp
    thermo_timint.cpp
    thermo_timint_expl.cpp
    thermo_timint_expleuler.cpp
    thermo_timint_genalpha.cpp
    thermo_timint_impl.cpp
    thermo_timint_ost.cpp
    thermo_timint_statics.cpp
    )

baci_add_library(alg_thermo ${OBJS_THERMO})

set(LIBS_ALG_THERMO
    # cmake-format: sortable
    adapter
    ele_thermo
    inpar
    )

baci_add_dependency(alg_thermo ${LIBS_ALG_THERMO})
