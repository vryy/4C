set(OBJS_DRT_RED_AIRWAYS_LIB
    # cmake-format: sortable
    acinus_evaluate.cpp
    acinus_impl.cpp
    airway_evaluate.cpp
    airway_impl.cpp
    airwayimplicitintegration.cpp
    inter_acinar_dep_evaluate.cpp
    inter_acinar_dep_impl.cpp
    red_acinus.cpp
    red_air_blood_scatra.cpp
    red_air_blood_scatra_evaluate.cpp
    red_air_blood_scatra_impl.cpp
    red_air_blood_scatraLine3.cpp
    red_air_blood_scatraLine3_evaluate.cpp
    red_air_blood_scatraLine3_impl.cpp
    red_airway.cpp
    red_airway_input.cpp
    red_airway_resulttest.cpp
    red_airways_dyn_drt.cpp
    red_inter_acinar_dep.cpp
    redairway_tissue.cpp
    )

baci_add_library(drt_red_airways ${OBJS_DRT_RED_AIRWAYS_LIB})

set(LIBS_DRT_RED_AIRWAYS
    # cmake-format: sortable
    drt_mat
    )

target_link_libraries(drt_red_airways PUBLIC ${LIBS_DRT_RED_AIRWAYS})
