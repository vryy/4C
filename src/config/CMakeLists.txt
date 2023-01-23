# configure a header file to pass some of the CMake settings
# to the source code
configure_file(
  "config_compile_settings.H.in" "${CMAKE_CURRENT_BINARY_DIR}/config_compile_settings.H"
  )

# get current Baci revision from git
include(GetGitRevision)
get_git_revision_information()
# configure revision file to pass git revision information to the source code
configure_file("config_revision.H.in" "${CMAKE_CURRENT_BINARY_DIR}/config_revision.H")

get_trilinos_version()
# configure trilinos version file to pass information to the source code
configure_file(
  "config_trilinos_version.H.in" "${CMAKE_CURRENT_BINARY_DIR}/config_trilinos_version.H"
  )

add_library(config INTERFACE)
target_include_directories(config INTERFACE "${CMAKE_CURRENT_BINARY_DIR}")
