# Exclude a file from unity build for its target.
function(four_c_unity_build_compile_separately _target _file)
  set_source_files_properties(
    ${_file} TARGET_DIRECTORY ${_target}_objs PROPERTIES SKIP_UNITY_BUILD_INCLUSION 1
    )
endfunction()
