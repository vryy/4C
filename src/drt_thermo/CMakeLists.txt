set(OBJS_DRT_THERMO_ELE
    # cmake-format: sortable
    thermo_ele_boundary_impl.cpp
    thermo_ele_impl.cpp
    thermo_element.cpp
    thermo_element_boundary_evaluate.cpp
    thermo_element_evaluate.cpp
    thermo_element_input.cpp
    )

baci_add_library(ele_thermo ${OBJS_DRT_THERMO_ELE})

set(LIBS_ELE_THERMO
    # cmake-format: sortable
    drt_inpar
    drt_lib
    )

target_link_libraries(ele_thermo PUBLIC ${LIBS_ELE_THERMO})

set(OBJS_DRT_THERMO
    # cmake-format: sortable
    thr_aux.cpp
    thr_dyn.cpp
    thr_resulttest.cpp
    thrtimint.cpp
    thrtimint_expl.cpp
    thrtimint_expleuler.cpp
    thrtimint_genalpha.cpp
    thrtimint_impl.cpp
    thrtimint_ost.cpp
    thrtimint_statics.cpp
    )

baci_add_library(alg_thermo ${OBJS_DRT_THERMO})

set(LIBS_ALG_THERMO
    # cmake-format: sortable
    drt_adapter
    drt_inpar
    ele_thermo
    )

target_link_libraries(alg_thermo PUBLIC ${LIBS_ALG_THERMO})
