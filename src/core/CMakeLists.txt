# Create the core module which makes the whole core functionality available.
baci_define_module(core)
add_library(baci::core ALIAS core)

target_link_libraries(core INTERFACE coupling discretization utils)
add_subdirectory(coupling)
add_subdirectory(discretization)
add_subdirectory(utils)
