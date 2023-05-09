set(OBJS_CREATERTD_MAIN
    # cmake-format: sortable
    create_rtdfiles_base.cpp
    )
set(OBJS_CREATERTD_LIB create_rtdfiles_wrapper.cpp create_rtdfiles_utils.cpp)
baci_add_library(create_rtdfiles ${OBJS_CREATERTD_LIB})
baci_add_dependency(create_rtdfiles lib)
# the create_rtdfiles executable
add_executable(create_rtd EXCLUDE_FROM_ALL ${OBJS_CREATERTD_MAIN})
baci_add_dependency(create_rtd create_rtdfiles)
set_target_properties(create_rtd PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_link_libraries(create_rtd PRIVATE -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group)

add_custom_target(rtd DEPENDS create_rtd)
