##
# Process global options
#

include(CheckCXXCompilerFlag)

#
# Add a compiler flag if the compiler understands it.
#
function(enable_flag_if_supported _flag)
  check_cxx_compiler_flag(${_flag} BACI_COMPILER_HAS_FLAG_${_flag})
  if(BACI_COMPILER_HAS_FLAG_${_flag})
    target_compile_options(baci_global_compile_settings INTERFACE ${_flag})
  endif()
endfunction()

# Define a target which pulls in all the global compiler settings and definitions
add_library(baci_global_compile_settings INTERFACE)

enable_flag_if_supported("-Wall")
enable_flag_if_supported("-Wextra")
# Disable unused parameter detection since there would be too many hits to fix
enable_flag_if_supported("-Wno-unused-parameter")
# Disable maybe uninitialized detection since this can give false positives that are hard to circumvent.
# Our address sanitizer checks will also find such errors, although only later.
enable_flag_if_supported("-Wno-maybe-uninitialized")
enable_flag_if_supported("-Wvla")

# For clang: do not error for a number of checks that are not yet fixed
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  enable_flag_if_supported("-Wno-error=overloaded-virtual")
  enable_flag_if_supported("-Wno-error=unused-variable")
  enable_flag_if_supported("-Wno-error=undefined-var-template")
  enable_flag_if_supported("-Wno-error=potentially-evaluated-expression")
endif()

baci_process_global_option(BACI_DSERROR_DUMP "dserror creates a core file" OFF)
baci_process_global_option(BACI_TRAP_FE "Crash BACI if a nan or inf occurs" ON)

baci_process_global_option(
  BACI_DEBUG
  "Turn on assertions and debug sections in code. Automatically turned on for DEBUG as CMAKE_BUILD_TYPE."
  OFF
  )
if(${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
  set(BACI_DEBUG "ON")
endif()
