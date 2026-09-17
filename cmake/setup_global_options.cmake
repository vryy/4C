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

# Backwards compatibility: if BUILD_SHARED_LIBS is set but not FOUR_C_BUILD_SHARED_LIBS, set the latter to the former.
if(NOT DEFINED FOUR_C_BUILD_SHARED_LIBS AND DEFINED BUILD_SHARED_LIBS)
  set(FOUR_C_BUILD_SHARED_LIBS
      ${BUILD_SHARED_LIBS}
      CACHE BOOL "Build shared libraries instead of static ones" FORCE
      )
  message(
    WARNING
      "Set FOUR_C_BUILD_SHARED_LIBS instead of BUILD_SHARED_LIBS. "
      "Setting FOUR_C_BUILD_SHARED_LIBS based on BUILD_SHARED_LIBS to ${FOUR_C_BUILD_SHARED_LIBS}."
    )
endif()

four_c_process_global_option(
  FOUR_C_BUILD_SHARED_LIBS
  DESCRIPTION
  "Build shared libraries instead of static ones"
  DEFAULT
  ON
  )
set(BUILD_SHARED_LIBS
    ${FOUR_C_BUILD_SHARED_LIBS}
    CACHE
      BOOL
      "Build shared libraries instead of static ones. Forces to be in sync with FOUR_C_BUILD_SHARED_LIBS."
      FORCE
    )

set(FOUR_C_VERIFY_INTERFACE_HEADER_SETS
    ${CMAKE_VERIFY_INTERFACE_HEADER_SETS}
    CACHE
      BOOL
      "Build in verify header mode. Forces to be in sync with CMAKE_VERIFY_INTERFACE_HEADER_SETS."
      FORCE
    )

four_c_process_global_option(
  FOUR_C_ENABLE_DEVELOPER_MODE
  DESCRIPTION
  "Enable developer mode (tries to optimize setup for iterative development cycles)"
  DEFAULT
  OFF
  )

if(MSVC)
  #if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang"
  #  AND NOT CMAKE_CXX_COMPILER_FRONTEND_VARIANT MATCHES "MSVC")
  # Enables the C++ standard-compliant preprocessor in MSVC
  #enable_compiler_flag_if_supported("/Zc:preprocessor")
  #endif()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    enable_compiler_flag_if_supported("/Zc:preprocessor")
  endif()

  # Enable high warning level (closest equivalent to -Wall -Wextra)
  enable_compiler_flag_if_supported("/W4")

  # Disable unused parameter warning (MSVC C4100)
  enable_compiler_flag_if_supported("/wd4100")

  # Disable overloaded virtual function warning (MSVC C4263 / C4264)
  enable_compiler_flag_if_supported("/wd4263")
  enable_compiler_flag_if_supported("/wd4264")

  # Note: Variable Length Arrays (-Wvla) are not supported by MSVC's C compiler,
  # and in C++ they are non-standard, so no flag is needed.

  # Note: -rdynamic is a Linux-specific linker flag. On Windows/MSVC,
  # visual studio exports symbols via CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS or module definition files.
else()
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
endif()

# Enable position-independent code. This flag is necessary to build shared libraries. Since our internal targets are not
# real libraries but object libraries, this property needs to be set explicitly.
if(FOUR_C_BUILD_SHARED_LIBS)
  message(VERBOSE "Enabling POSITION_INDEPENDENT_CODE on internal targets.")
  set_target_properties(
    four_c_private_compile_interface PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE TRUE
    )
endif()

# Special flags for GCC
message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # Disable maybe uninitialized detection since this can give false positives that are hard to circumvent.
  # Our address sanitizer checks will also find such errors, although only later.
  enable_compiler_flag_if_supported("-Wno-maybe-uninitialized")
endif()

four_c_process_global_option(
  FOUR_C_ENABLE_WARNINGS_AS_ERRORS
  DESCRIPTION
  "Treat warnings as errors when compiling"
  DEFAULT
  OFF
  )
if(FOUR_C_ENABLE_WARNINGS_AS_ERRORS)
  if(MSVC)
    enable_compiler_flag_if_supported("/WX")
  else()
    enable_compiler_flag_if_supported("-Werror")
  endif()
endif()

four_c_process_global_option(
  FOUR_C_ENABLE_NATIVE_OPTIMIZATIONS
  DESCRIPTION
  "Optimize for current hardware"
  DEFAULT
  OFF
  )
if(FOUR_C_ENABLE_NATIVE_OPTIMIZATIONS)
  enable_compiler_flag_if_supported("-march=native")
endif()

# If enabled, build all targets with address sanitizer
four_c_process_global_option(
  FOUR_C_ENABLE_ADDRESS_SANITIZER
  DESCRIPTION
  "Compile with address sanitizer"
  DEFAULT
  OFF
  )
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
  FOUR_C_ENABLE_COVERAGE
  DESCRIPTION
  "Set up a build to gather coverage information with LLVM source based coverage"
  DEFAULT
  OFF
  )
if(FOUR_C_ENABLE_COVERAGE)
  four_c_check_compiles(
    FOUR_C_COMPILER_SUPPORT_COVERAGE
    COMPILE_OPTIONS
    "-fprofile-instr-generate"
    "-fcoverage-mapping"
    LINK_OPTIONS
    "-fprofile-instr-generate"
    "-fcoverage-mapping"
    "-Wl,--build-id=sha1"
    APPEND_ON_SUCCESS
    )

  if(NOT FOUR_C_COMPILER_SUPPORT_COVERAGE)
    message(
      FATAL_ERROR
        "Option FOUR_C_ENABLE_COVERAGE is ON but the compiler does not support this feature."
      )
  endif()
endif()

four_c_process_global_option(
  FOUR_C_ENABLE_CORE_DUMP
  DESCRIPTION
  "Uncaught exceptions create a core file"
  DEFAULT
  OFF
  )

four_c_process_global_option(
  FOUR_C_ENABLE_FE_TRAPPING
  DESCRIPTION
  "Crash the program if a nan or inf would occur"
  DEFAULT
  ON
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
  if(MSVC)
    enable_compiler_flag_if_supported("/fp:fast")
  else()
    enable_compiler_flag_if_supported("-fno-trapping-math")
  endif()
endif()

four_c_process_global_option(
  FOUR_C_ENABLE_IWYU
  DESCRIPTION
  "Enable include-what-you-use"
  DEFAULT
  OFF
  )
if(FOUR_C_ENABLE_IWYU)
  find_program(FOUR_C_IWYU_EXECUTABLE NAMES include-what-you-use iwyu)
  if(NOT FOUR_C_IWYU_EXECUTABLE)
    message(
      FATAL_ERROR "Option FOUR_C_ENABLE_IWYU is ON but include-what-you-use/iwyu is not found."
                  "You can specify the path in the CMake variable FOUR_C_IWYU_EXECUTABLE."
      )
  endif()
endif()

four_c_process_global_option(
  FOUR_C_ENABLE_PYTHON_BINDINGS
  DESCRIPTION
  "Enable building the 4C Python bindings (py4C) using pybind11."
  DEFAULT
  OFF
  )

four_c_process_global_option(
  FOUR_C_CLANGCUDA
  DESCRIPTION
  "Enable the relevant CMake compile definitions needed to use utilities/clangcuda++ as the compiler. This is currently necessary to use the CUDA backend of Kokkos in 4C, e.g. along with MIRCO."
  DEFAULT
  OFF
  )
if(FOUR_C_CLANGCUDA)
  if(FOUR_C_WITH_ARBORX)
    message(
      WARNING
        "Enabling both FOUR_C_CLANGCUDA and FOUR_C_WITH_ARBORX is not advised. This requires using an external CUDA-enabled ArborX installation and has not been tested."
      )
  endif()
  if(FOUR_C_WITH_VTK)
    message(
      FATAL_ERROR
        "Enabling both FOUR_C_CLANGCUDA and FOUR_C_WITH_VTK currently causes compilation errors."
      )
  endif()
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

# using /Z7 instead of /Zi because /Zi prevents ccache to work correctly
set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "Embedded")

if(${FOUR_C_BUILD_TYPE_UPPER} MATCHES DEBUG)
  set(FOUR_C_ENABLE_ASSERTIONS
      "ON"
      CACHE BOOL "Forced ON due to build type DEBUG" FORCE
      )
  if(MSVC)
    target_compile_options(four_c_private_compile_interface INTERFACE "/Od")
    target_link_options(four_c_private_compile_interface INTERFACE "/DEBUG")
  else()
    target_compile_options(four_c_private_compile_interface INTERFACE "-O0")
    target_link_options(four_c_private_compile_interface INTERFACE "-O0")

    target_compile_options(four_c_private_compile_interface INTERFACE "-g")
  endif()
endif()

if(${FOUR_C_BUILD_TYPE_UPPER} MATCHES RELEASE)
  if(MSVC)
    target_compile_options(four_c_private_compile_interface INTERFACE "/O2" "/Ob2")
  else()
    target_compile_options(four_c_private_compile_interface INTERFACE "-O3")
    target_link_options(four_c_private_compile_interface INTERFACE "-O3")

    enable_compiler_flag_if_supported("-funroll-loops")
  endif()
endif()

if(${FOUR_C_BUILD_TYPE_UPPER} MATCHES RELWITHDEBINFO)
  if(MSVC)
    target_compile_options(four_c_private_compile_interface INTERFACE "/O2" "/Ob1")
    target_link_options(four_c_private_compile_interface INTERFACE "/DEBUG" "/LTCG")
  else()
    target_compile_options(four_c_private_compile_interface INTERFACE "-O3")
    target_link_options(four_c_private_compile_interface INTERFACE "-O3")

    target_compile_options(four_c_private_compile_interface INTERFACE "-g")
    enable_compiler_flag_if_supported("-funroll-loops")
  endif()
endif()

# Evaluate this option now to get the correct output in case it is forced ON in DEBUG mode.
four_c_process_global_option(
  FOUR_C_ENABLE_ASSERTIONS
  DESCRIPTION
  "Turn on assertions and debug sections in 4C code, and turn on assertions in the standard library. Automatically turned on for DEBUG as CMAKE_BUILD_TYPE."
  DEFAULT
  OFF
  )

if(FOUR_C_ENABLE_ASSERTIONS)
  # In case we enable assertions, ensure that we also enable the standard library assertions.
  # Note that starting from GCC 15, they are turned on in unoptimized -O0 builds. Since we
  # usually use optimization in testing, we need to explicitly turn them on when requesting
  # 4C's own assertions.
  # Note that we can happily define this libstdc++ macro regardless of the standard library
  # implementation, as it will just be ignored by other implementations.
  target_compile_options(four_c_private_compile_interface INTERFACE "-D_GLIBCXX_ASSERTIONS")
endif()

four_c_process_global_option(
  FOUR_C_ENABLE_METADATA_GENERATION
  DESCRIPTION
  "Generate metadata after building 4C. Requires python and invokes 4C."
  DEFAULT
  ON
  )

##
# Add potential user flags at the very end
##

four_c_process_cache_variable(
  FOUR_C_CXX_FLAGS
  TYPE
  STRING
  DESCRIPTION
  "Expert setting. Additional C++ compiler flags to use for all 4C targets. Note that these flags are added at the very end and can override various other flags."
  DEFAULT
  ""
  )
four_c_process_cache_variable(
  FOUR_C_CXX_LINKER_FLAGS
  TYPE
  STRING
  DESCRIPTION
  "Expert setting. Additional C++ linker flags to use for all 4C targets. Note that these flags are added at the very end and can override various other flags."
  DEFAULT
  ""
  )

# Compiler
separate_arguments(_split UNIX_COMMAND ${FOUR_C_CXX_FLAGS})
target_compile_options(four_c_private_compile_interface INTERFACE ${_split})

# Linker
separate_arguments(_split UNIX_COMMAND ${FOUR_C_CXX_LINKER_FLAGS})
target_link_options(four_c_private_compile_interface INTERFACE ${_split})

### Do not add any more flags here! User flags have already been added.

message(STATUS "") # Separate output with empty line
