set(OBJS_RED_AIRWAYS_LIB
    # cmake-format: sortable
    red_airways_acinus.cpp
    red_airways_acinus_evaluate.cpp
    red_airways_acinus_impl.cpp
    red_airways_air_blood_scatra.cpp
    red_airways_air_blood_scatra_evaluate.cpp
    red_airways_air_blood_scatra_impl.cpp
    red_airways_air_blood_scatraLine3.cpp
    red_airways_air_blood_scatraLine3_evaluate.cpp
    red_airways_air_blood_scatraLine3_impl.cpp
    red_airways_airway.cpp
    red_airways_airway_evaluate.cpp
    red_airways_airway_impl.cpp
    red_airways_dyn_drt.cpp
    red_airways_implicitintegration.cpp
    red_airways_input.cpp
    red_airways_interacinardep.cpp
    red_airways_interacinardep_evaluate.cpp
    red_airways_interacinardep_impl.cpp
    red_airways_resulttest.cpp
    red_airways_tissue.cpp
    )

baci_add_library(red_airways ${OBJS_RED_AIRWAYS_LIB})

set(LIBS_RED_AIRWAYS
    # cmake-format: sortable
    mat
    )

baci_add_dependency(red_airways ${LIBS_RED_AIRWAYS})
