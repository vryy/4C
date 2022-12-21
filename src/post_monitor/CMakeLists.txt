set(OBJS_POST_MONITOR
    # cmake-format: sortable
    post_monitor.cpp
    )

add_executable(post_monitor EXCLUDE_FROM_ALL ${OBJS_POST_MONITOR})
set_target_properties(post_monitor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_link_libraries(post_monitor PUBLIC ${FILTER_LIBRARIES})
