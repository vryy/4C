# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

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
  four_c_check_compiles(
    FOUR_C_COMPILER_HAS_FLAG_${_flag_var} COMPILE_OPTIONS ${_flag} APPEND_ON_SUCCESS
    )
endfunction()

#
# Add a linker flag if the linker understands it.
#
function(enable_linker_flag_if_supported _flag)
  # Clean up a flag to get a clean variable name. The flag is passed to the
  # compiler as a define.
  string(REGEX REPLACE "^-" "" _flag_var ${_flag})
  string(REGEX REPLACE "\[-=\]" "_" _flag_var ${_flag_var})
  four_c_check_compiles(FOUR_C_LINKER_HAS_FLAG_${_flag_var} LINK_OPTIONS ${_flag} APPEND_ON_SUCCESS)
endfunction()

# Enable all warnings that are supported by the compiler
enable_compiler_flag_if_supported("-Wall")
enable_compiler_flag_if_supported("-Wextra")
enable_compiler_flag_if_supported("-Wvla")

# Disable unused parameter detection since there would be too many hits to fix
enable_compiler_flag_if_supported("-Wno-unused-parameter")
# Disable overloaded virtual function detection. This requires a lot of architectural changes to fix.
enable_compiler_flag_if_supported("-Wno-overloaded-virtual")

# Export symbols (necessary for stacktraces)
enable_linker_flag_if_supported("-rdynamic")

# Enable position-independent code. This flag is necessary to build shared libraries. Since our internal targets are not
# real libaries but object libraries, this property needs to be set explicitly.
if(BUILD_SHARED_LIBS)
  message(VERBOSE "Enabling POSITION_INDEPENDENT_CODE on internal targets.")
  set_target_properties(
    four_c_private_compile_interface PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE TRUE
    )
endif()

# Special flags for GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # Disable maybe uninitialized detection since this can give false positives that are hard to circumvent.
  # Our address sanitizer checks will also find such errors, although only later.
  enable_compiler_flag_if_supported("-Wno-maybe-uninitialized")
endif()

# For clang: do not error for a number of checks that are not yet fixed
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  enable_compiler_flag_if_supported("-Wno-error=templates")
endif()

four_c_process_global_option(
  FOUR_C_ENABLE_WARNINGS_AS_ERRORS "Treat warnings as errors when compiling" OFF
  )
if(FOUR_C_ENABLE_WARNINGS_AS_ERRORS)
  enable_compiler_flag_if_supported("-Werror")
endif()

four_c_process_global_option(FOUR_C_ENABLE_NATIVE_OPTIMIZATIONS "Optimize for current hardware" OFF)
if(FOUR_C_ENABLE_NATIVE_OPTIMIZATIONS)
  enable_compiler_flag_if_supported("-march=native")
endif()

# If enabled, build all targets with address sanitizer
four_c_process_global_option(FOUR_C_ENABLE_ADDRESS_SANITIZER "Compile with address sanitizer" OFF)
if(FOUR_C_ENABLE_ADDRESS_SANITIZER)
  # We get better stack traces in ASAN with this flag.
  enable_compiler_flag_if_supported("-fno-omit-frame-pointer")
  enable_linker_flag_if_supported("-fno-omit-frame-pointer")

  # ASAN requires to check both flags at once.
  four_c_check_compiles(
    FOUR_C_COMPILER_LINKER_SUPPORT_ASAN
    COMPILE_OPTIONS
    "-fsanitize=address"
    LINK_OPTIONS
    "-fsanitize=address"
    APPEND_ON_SUCCESS
    )

  if(NOT FOUR_C_COMPILER_LINKER_SUPPORT_ASAN)
    message(
      FATAL_ERROR
        "Option FOUR_C_ENABLE_ADDRESS_SANITIZER is ON but the compiler does not support this feature."
      )
  endif()
endif()

four_c_process_global_option(
  FOUR_C_ENABLE_COVERAGE "Set up a build to gather coverage information" OFF
  )
if(FOUR_C_ENABLE_COVERAGE)
  four_c_check_compiles(
    FOUR_C_COMPILER_SUPPORT_COVERAGE
    COMPILE_OPTIONS
    "-fprofile-arcs"
    "-ftest-coverage"
    LINK_OPTIONS
    "-fprofile-arcs"
    "-ftest-coverage"
    APPEND_ON_SUCCESS
    )

  if(NOT FOUR_C_COMPILER_SUPPORT_COVERAGE)
    message(
      FATAL_ERROR
        "Option FOUR_C_ENABLE_COVERAGE is ON but the compiler does not support this feature."
      )
  endif()
endif()

four_c_process_global_option(FOUR_C_ENABLE_CORE_DUMP "Uncaught exceptions create a core file" OFF)

four_c_process_global_option(
  FOUR_C_ENABLE_FE_TRAPPING "Crash the program if a nan or inf would occur" ON
  )
# We need to let the compiler know that we intend to use the floating-point trapping mechanism.
if(FOUR_C_ENABLE_FE_TRAPPING)
  enable_compiler_flag_if_supported("-ftrapping-math")
  # Let's just make this an error to avoid very hard to debug errors.
  if(NOT FOUR_C_COMPILER_HAS_FLAG_ftrapping_math)
    message(
      FATAL_ERROR
        "Option FOUR_C_ENABLE_FE_TRAPPING is ON but the compiler does not support this feature."
        "Specifically, the compiler does not support -ftrapping-math, which is necessary to"
        "generate code that can safely use the floating-point trapping mechanism."
      )
  endif()
else()
  enable_compiler_flag_if_supported("-fno-trapping-math")
endif()

##
# Optimization flags
# These flags are reasonable defaults. Users may amend them by setting FOUR_C_CXX_FLAGS and/or FOUR_C_CXX_FLAGS_<CONFIG>.
# The default way to add more flags via CMAKE_CXX_FLAGS(_<CONFIG>) is also supported but note that these flags are
# added in front of all other flags and cannot override the defaults set below.
#
# For build types DEBUG, RELEASE, and RELWITHDEBINFO, Cmake populates the flags with default optimization levels.
# For our own code, we explicitly set the flags again, to be very clear in what we use for compilation. Projects included
# via fetch_content() will just use the defaults plus whatever is added to CMAKE_CXX_FLAGS(_<CONFIG>).
#
# Side note: we use the same flag for a few legacy C sources, which we compile as C++ anyway.
##

if(${FOUR_C_BUILD_TYPE_UPPER} MATCHES DEBUG)
  set(FOUR_C_ENABLE_ASSERTIONS
      "ON"
      CACHE BOOL "Forced ON due to build type DEBUG" FORCE
      )
  target_compile_options(four_c_private_compile_interface INTERFACE "-O0")
  target_link_options(four_c_private_compile_interface INTERFACE "-O0")

  target_compile_options(four_c_private_compile_interface INTERFACE "-g")
endif()

if(${FOUR_C_BUILD_TYPE_UPPER} MATCHES RELEASE)
  target_compile_options(four_c_private_compile_interface INTERFACE "-O3")
  target_link_options(four_c_private_compile_interface INTERFACE "-O3")

  enable_compiler_flag_if_supported("-funroll-loops")
endif()

if(${FOUR_C_BUILD_TYPE_UPPER} MATCHES RELWITHDEBINFO)
  target_compile_options(four_c_private_compile_interface INTERFACE "-O2")
  target_link_options(four_c_private_compile_interface INTERFACE "-O2")

  target_compile_options(four_c_private_compile_interface INTERFACE "-g")
  enable_compiler_flag_if_supported("-funroll-loops")
endif()

# Evaluate this option now to get the correct output in case it is force ON in DEBUG mode.
four_c_process_global_option(
  FOUR_C_ENABLE_ASSERTIONS
  "Turn on assertions and debug sections in code. Automatically turned on for DEBUG as CMAKE_BUILD_TYPE."
  OFF
  )

##
# Add potential user flags at the very end
##

# Compiler
separate_arguments(_split UNIX_COMMAND ${FOUR_C_CXX_FLAGS})
target_compile_options(four_c_private_compile_interface INTERFACE ${_split})
separate_arguments(_split UNIX_COMMAND ${FOUR_C_CXX_FLAGS_${FOUR_C_BUILD_TYPE_UPPER}})
target_compile_options(four_c_private_compile_interface INTERFACE ${_split})

# Linker
separate_arguments(_split UNIX_COMMAND ${FOUR_C_CXX_LINKER_FLAGS})
target_link_options(four_c_private_compile_interface INTERFACE ${_split})
separate_arguments(_split UNIX_COMMAND ${FOUR_C_CXX_LINKER_FLAGS_${FOUR_C_BUILD_TYPE_UPPER}})
target_link_options(four_c_private_compile_interface INTERFACE ${_split})

### Do not add any more flags here! User flags have already been added.
