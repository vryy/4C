set(OBJS_POROMULTIPHASE_SCATRA
    # cmake-format: sortable
    poromultiphase_scatra_base.cpp
    poromultiphase_scatra_dyn.cpp
    poromultiphase_scatra_function.cpp
    poromultiphase_scatra_monolithic_twoway.cpp
    poromultiphase_scatra_partitioned_twoway.cpp
    poromultiphase_scatra_utils.cpp
    )

baci_add_library(alg_poromultiphase_scatra ${OBJS_POROMULTIPHASE_SCATRA})

set(LIBS_ALG_PROMULTIPHASE_SCATRA
    # cmake-format: sortable
    alg_poromultiphase
    )

target_link_libraries(alg_poromultiphase_scatra PUBLIC ${LIBS_ALG_PROMULTIPHASE_SCATRA})

set(OBJS_POROMULTIPHASE_ART_COUPL
    # cmake-format: sortable
    poromultiphase_scatra_artery_coupling_base.cpp
    poromultiphase_scatra_artery_coupling_linebased.cpp
    poromultiphase_scatra_artery_coupling_nodebased.cpp
    poromultiphase_scatra_artery_coupling_nodetopoint.cpp
    poromultiphase_scatra_artery_coupling_nonconforming.cpp
    poromultiphase_scatra_artery_coupling_pair.cpp
    poromultiphase_scatra_artery_coupling_surfbased.cpp
    )

baci_add_library(lib_poromulti_art_coupl ${OBJS_POROMULTIPHASE_ART_COUPL})

set(LIBS_LIB_POROMULTI_ART_COUPL
    # cmake-format: sortable
    alg_artnet
    linalg
    )

target_link_libraries(lib_poromulti_art_coupl PUBLIC ${LIBS_LIB_POROMULTI_ART_COUPL})
