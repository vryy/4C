set(OBJS_DRT_THERMO_ELE
    thermo_ele_boundary_impl.cpp
    thermo_ele_impl.cpp
    thermo_element.cpp
    thermo_element_boundary_evaluate.cpp
    thermo_element_evaluate.cpp
    thermo_element_input.cpp
    )

baci_add_library(ele_thermo ${OBJS_DRT_THERMO_ELE})
target_link_libraries(ele_thermo PUBLIC drt_lib drt_inpar)

set(OBJS_DRT_THERMO
    thr_aux.cpp
    thr_dyn.cpp
    thr_resulttest.cpp
    thrtimint.cpp
    thrtimint_genalpha.cpp
    thrtimint_impl.cpp
    thrtimint_ost.cpp
    thrtimint_statics.cpp
    thrtimint_expl.cpp
    thrtimint_expleuler.cpp
    )

baci_add_library(alg_thermo ${OBJS_DRT_THERMO})
target_link_libraries(alg_thermo PUBLIC drt_inpar drt_adapter)
