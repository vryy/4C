set(OBJS_ARTNET
    # cmake-format: sortable
    art_net_dyn_drt.cpp
    art_net_impl_stationary.cpp
    art_net_timint.cpp
    art_net_utils.cpp
    artnetexplicitintegration.cpp
    )

baci_add_library(alg_artnet ${OBJS_ARTNET})

set(LIBS_ALG_ARTNET
    # cmake-format: sortable
    adapter
    adapter
    ele_scatra
    )
target_link_libraries(alg_artnet PUBLIC ${LIBS_ALG_ARTNET})

set(OBJS_ARTNET_LIB
    # cmake-format: sortable
    art_junction.cpp
    art_terminal_bc.cpp
    art_write_gnuplot.cpp
    artery.cpp
    artery_ele_calc.cpp
    artery_ele_calc_lin_exp.cpp
    artery_ele_calc_pres_based.cpp
    artery_ele_factory.cpp
    artery_evaluate.cpp
    artery_input.cpp
    artery_resulttest.cpp
    )

baci_add_library(lib_artnet ${OBJS_ARTNET_LIB})

set(LIBS_LIB_ARTNET
    # cmake-format: sortable
    adapter
    ele_scatra
    inpar
    lib
    mat
    )

target_link_libraries(lib_artnet PUBLIC ${LIBS_LIB_ARTNET})
