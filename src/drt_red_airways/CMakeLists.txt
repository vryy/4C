set(OBJS_DRT_RED_AIRWAYS_LIB
    acinus_evaluate.cpp
    airway_evaluate.cpp
    inter_acinar_dep_evaluate.cpp
    red_air_blood_scatra_evaluate.cpp
    red_air_blood_scatraLine3_evaluate.cpp
    airway_impl.cpp
    acinus_impl.cpp
    inter_acinar_dep_impl.cpp
    red_air_blood_scatra_impl.cpp
    red_air_blood_scatraLine3_impl.cpp
    red_airway.cpp
    red_acinus.cpp
    red_inter_acinar_dep.cpp
    red_air_blood_scatra.cpp
    red_air_blood_scatraLine3.cpp
    red_airway_input.cpp
    red_airway_resulttest.cpp
    airwayimplicitintegration.cpp
    redairway_tissue.cpp
    red_airways_dyn_drt.cpp
    )

baci_add_library(drt_red_airways ${OBJS_DRT_RED_AIRWAYS_LIB})
target_link_libraries(drt_red_airways PUBLIC drt_mat)
