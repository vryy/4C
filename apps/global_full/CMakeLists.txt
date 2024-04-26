set(OBJS_FOUR_C_MAIN
    # cmake-format: sortable
    4C_global_full_cal_control.cpp
    4C_global_full_control.cpp
    4C_global_full_init_control.cpp
    4C_global_full_inp_control.cpp
    4C_global_full_main.cpp
    )

add_executable(${FOUR_C_EXECUTABLE_NAME} ${OBJS_FOUR_C_MAIN})
set_target_properties(
  ${FOUR_C_EXECUTABLE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
four_c_set_up_executable(${FOUR_C_EXECUTABLE_NAME})

add_custom_target(baci DEPENDS ${FOUR_C_EXECUTABLE_NAME})
