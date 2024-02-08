# configure a header file to pass some of the CMake settings
# to the source code
configure_file("baci_config.H.in" "${CMAKE_CURRENT_BINARY_DIR}/baci_config.H")

baci_get_git_revision_information()
# configure revision file to pass git revision information to the source code
configure_file("baci_config_revision.H.in" "${CMAKE_CURRENT_BINARY_DIR}/baci_config_revision.H")

baci_get_trilinos_version()
# configure trilinos version file to pass information to the source code
configure_file(
  "baci_config_trilinos_version.H.in" "${CMAKE_CURRENT_BINARY_DIR}/baci_config_trilinos_version.H"
  )

set(_headers
    ${CMAKE_CURRENT_BINARY_DIR}/baci_config.H
    ${CMAKE_CURRENT_BINARY_DIR}/baci_config_revision.H
    ${CMAKE_CURRENT_BINARY_DIR}/baci_config_trilinos_version.H
    )

baci_auto_define_module()

# Manually add the configured headers to the config target
target_sources(
  ${AUTO_DEFINED_MODULE_NAME}_deps
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_BINARY_DIR}
            FILES
            ${_headers}
  )
