set(OBJS_CREATERTD
    4C_create_rtdfiles_base.cpp 4C_create_rtdfiles_wrapper.cpp 4C_create_rtdfiles_utils.cpp
    )
add_executable(create_rtd EXCLUDE_FROM_ALL ${OBJS_CREATERTD})
set_target_properties(create_rtd PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
four_c_set_up_executable(create_rtd)
add_custom_target(rtd DEPENDS create_rtd)
