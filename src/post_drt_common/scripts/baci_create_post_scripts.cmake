function(copy_script script_name)
  add_custom_command(
    TARGET post_processor
    POST_BUILD
    COMMAND
      ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${script_name}
      ${PROJECT_BINARY_DIR}
    COMMENT "Create script ${script_name}"
    )
endfunction()

copy_script(post_drt_ensight)
copy_script(post_drt_gid)
copy_script(post_drt_vti)
copy_script(post_drt_vtu)
