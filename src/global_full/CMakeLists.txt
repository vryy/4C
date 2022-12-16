set(OBJS_BACI_MAIN
    # cmake-format: sortable
    baci.cpp
    global_cal_control.cpp
    global_control.cpp
    global_init_control.cpp
    global_inp_control.cpp
    )

# the baci executable
add_executable(${baciname} ${OBJS_BACI_MAIN})
set_target_properties(${baciname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_link_libraries(${baciname} ${BACI_LIBRARIES})

add_custom_target(baci DEPENDS ${baciname})
