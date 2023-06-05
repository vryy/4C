set(OBJS_BACI_MAIN
    # cmake-format: sortable
    global_full_baci.cpp
    global_full_cal_control.cpp
    global_full_control.cpp
    global_full_init_control.cpp
    global_full_inp_control.cpp
    )

# the baci executable
add_executable(${baciname} ${OBJS_BACI_MAIN})
set_target_properties(${baciname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_target_requires_cyclic_legacy_linking(${baciname})

add_custom_target(baci DEPENDS ${baciname})
