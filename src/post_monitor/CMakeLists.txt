set(OBJS_POST_DRT_MONITOR
    # cmake-format: sortable
    post_drt_monitor.cpp
    )

add_executable(post_drt_monitor EXCLUDE_FROM_ALL ${OBJS_POST_DRT_MONITOR})
set_target_properties(post_drt_monitor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_link_libraries(post_drt_monitor PUBLIC ${FILTER_LIBRARIES})
