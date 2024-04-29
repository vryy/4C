function(four_c_unity_build_compile_separately _file)
  set_source_files_properties(${_file} PROPERTIES SKIP_UNITY_BUILD_INCLUSION 1)
endfunction()
