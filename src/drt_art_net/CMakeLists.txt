set(OBJS_DRT_ARTNET
    art_net_dyn_drt.cpp
    art_net_utils.cpp
    art_net_timint.cpp
    art_net_impl_stationary.cpp
    artnetexplicitintegration.cpp
    )

baci_add_library(alg_artnet ${OBJS_DRT_ARTNET})
target_link_libraries(alg_artnet PUBLIC drt_adapter)

set(OBJS_DRT_ARTNET_LIB
    artery.cpp
    artery_ele_calc_lin_exp.cpp
    artery_ele_calc_pres_based.cpp
    artery_ele_factory.cpp
    artery_evaluate.cpp
    artery_ele_calc.cpp
    art_junction.cpp
    art_terminal_bc.cpp
    art_write_gnuplot.cpp
    artery_input.cpp
    artery_resulttest.cpp
    )

baci_add_library(drt_lib_artnet ${OBJS_DRT_ARTNET_LIB})
target_link_libraries(drt_lib_artnet PUBLIC drt_lib_inpar)
