##
# Process global options
#

#
# Add a compiler flag if the compiler understands it.
#
function(enable_compiler_flag_if_supported _flag)
  # Clean up a flag to get a clean variable name. The flag is passed to the
  # compiler as a define.
  string(REGEX REPLACE "^-" "" _flag_var ${_flag})
  string(REGEX REPLACE "\[-=\]" "_" _flag_var ${_flag_var})
  baci_add_settings_if_compiles(BACI_COMPILER_HAS_FLAG_${_flag_var} COMPILE_OPTIONS ${_flag})
endfunction()

#
# Add a linker flag if the linker understands it.
#
function(enable_linker_flag_if_supported _flag)
  # Clean up a flag to get a clean variable name. The flag is passed to the
  # compiler as a define.
  string(REGEX REPLACE "^-" "" _flag_var ${_flag})
  string(REGEX REPLACE "\[-=\]" "_" _flag_var ${_flag_var})
  baci_add_settings_if_compiles(BACI_LINKER_HAS_FLAG_${_flag_var} LINK_OPTIONS ${_flag})
endfunction()

enable_compiler_flag_if_supported("-Wall")
enable_compiler_flag_if_supported("-Wextra")
# Disable unused parameter detection since there would be too many hits to fix
enable_compiler_flag_if_supported("-Wno-unused-parameter")
# Disable maybe uninitialized detection since this can give false positives that are hard to circumvent.
# Our address sanitizer checks will also find such errors, although only later.
enable_compiler_flag_if_supported("-Wno-maybe-uninitialized")
enable_compiler_flag_if_supported("-Wvla")

# Export symbols (necessary for stacktraces)
enable_linker_flag_if_supported("-rdynamic")

# Enable position-independent code. This flag is necessary to build shared libraries. Since our internal targets are not
# real libaries but object libraries, this property needs to be set explicitly.
if(BUILD_SHARED_LIBS)
  message(VERBOSE "Enabling POSITION_INDEPENDENT_CODE on internal targets.")
  set_target_properties(
    baci_private_compile_interface PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE TRUE
    )
endif()

# For clang: do not error for a number of checks that are not yet fixed
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  enable_compiler_flag_if_supported("-Wno-error=overloaded-virtual")
  enable_compiler_flag_if_supported("-Wno-error=unused-variable")
  enable_compiler_flag_if_supported("-Wno-error=undefined-var-template")
  enable_compiler_flag_if_supported("-Wno-error=potentially-evaluated-expression")
endif()

# If enabled, build all targets with address sanitizer
option(BACI_WITH_ADDRESS_SANITIZER "Compile BACI with address sanitizer" OFF)
if(BACI_WITH_ADDRESS_SANITIZER)
  # We get better stack traces in ASAN with this flag.
  enable_compiler_flag_if_supported("-fno-omit-frame-pointer")
  enable_linker_flag_if_supported("-fno-omit-frame-pointer")

  # ASAN requires to check both flags at once.
  baci_add_settings_if_compiles(
    BACI_COMPILER_LINKER_SUPPORT_ASAN
    COMPILE_OPTIONS
    "-fsanitize=address"
    LINK_OPTIONS
    "-fsanitize=address"
    )

  if(NOT BACI_COMPILER_LINKER_SUPPORT_ASAN)
    message(
      FATAL_ERROR
        "Option BACI_WITH_ADDRESS_SANITIZER is ON but the compiler does not support this feature."
      )
  endif()
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
