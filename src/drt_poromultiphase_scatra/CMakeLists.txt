set(OBJS_DRT_POROMULTIPHASE_SCATRA
    poromultiphase_scatra_dyn.cpp
    poromultiphase_scatra_utils.cpp
    poromultiphase_scatra_base.cpp
    poromultiphase_scatra_partitioned_twoway.cpp
    poromultiphase_scatra_monolithic_twoway.cpp
    poromultiphase_scatra_function.cpp
    )

baci_add_library(alg_poromultiphase_scatra ${OBJS_DRT_POROMULTIPHASE_SCATRA})

target_link_libraries(alg_poromultiphase_scatra PUBLIC alg_poromultiphase)

set(OBJS_DRT_POROMULTIPHASE_ART_COUPL
    poromultiphase_scatra_artery_coupling_base.cpp
    poromultiphase_scatra_artery_coupling_nonconforming.cpp
    poromultiphase_scatra_artery_coupling_nodebased.cpp
    poromultiphase_scatra_artery_coupling_nodetopoint.cpp
    poromultiphase_scatra_artery_coupling_linebased.cpp
    poromultiphase_scatra_artery_coupling_surfbased.cpp
    poromultiphase_scatra_artery_coupling_pair.cpp
    )

baci_add_library(drt_lib_poromulti_art_coupl ${OBJS_DRT_POROMULTIPHASE_ART_COUPL})
target_link_libraries(drt_lib_poromulti_art_coupl PUBLIC linalg)
