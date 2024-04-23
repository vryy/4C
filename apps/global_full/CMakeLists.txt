set(OBJS_FOUR_C_MAIN
    # cmake-format: sortable
    baci_global_full_baci.cpp
    baci_global_full_cal_control.cpp
    baci_global_full_control.cpp
    baci_global_full_init_control.cpp
    baci_global_full_inp_control.cpp
    )

# the baci executable
add_executable(${baciname} ${OBJS_FOUR_C_MAIN})
set_target_properties(${baciname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_set_up_executable(${baciname})

add_custom_target(baci DEPENDS ${baciname})
