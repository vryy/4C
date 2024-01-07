##
# Process global options
#

# Define a target which pulls in all the global compiler settings and definitions
add_library(baci_global_compile_settings INTERFACE)

target_compile_options(baci_global_compile_settings INTERFACE "-Wall")
target_compile_options(baci_global_compile_settings INTERFACE "-Wextra")
# Disable unused parameter detection since there would be too many hits to fix
target_compile_options(baci_global_compile_settings INTERFACE "-Wno-unused-parameter")
target_compile_options(baci_global_compile_settings INTERFACE "-Wvla")

baci_process_global_option(DSERROR_DUMP "dserror creates a core file" OFF)
baci_process_global_option(TRAP_FE "Crash BACI if a nan or inf occurs" ON)

baci_process_global_option(
  BACI_DEBUG
  "Turn on assertions and debug sections in code. Automatically turned on for DEBUG as CMAKE_BUILD_TYPE."
  OFF
  )
if(${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
  set(BACI_DEBUG "ON")
endif()
