# configure a header file to pass some of the CMake settings
# to the source code
configure_file("4C_config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/4C_config.hpp")

baci_get_git_revision_information()
# configure revision file to pass git revision information to the source code
configure_file("4C_config_revision.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/4C_config_revision.hpp")

baci_get_trilinos_version()
# configure trilinos version file to pass information to the source code
configure_file(
  "4C_config_trilinos_version.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/4C_config_trilinos_version.hpp"
  )

set(_headers
    ${CMAKE_CURRENT_BINARY_DIR}/4C_config.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/4C_config_revision.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/4C_config_trilinos_version.hpp
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
