set(OBJS_POST_DRT_S8CONVERT
    # cmake-format: sortable
    post_drt_s8convert.cpp
    )

add_executable(post_drt_s8convert EXCLUDE_FROM_ALL ${OBJS_POST_DRT_S8CONVERT})
set_target_properties(post_drt_s8convert PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_link_libraries(post_drt_s8convert ${FILTER_LIBRARIES})
