#
# Find the backtrace includes and library
# backtrace is crucial for boost-stacktrace to print out the line numbers. This is documented here:
#   https://github.com/boostorg/stacktrace/issues/97
#
# There is a FindBacktrace.cmake in cmake default modules but it is useless/broken, hence we write our own module.
#
# Backtrace_INCLUDE_DIR - where to find backtrace.h
# Backtrace_LIBRARIES   - List of fully qualified libraries to link against.
# Backtrace_FOUND       - Do not attempt to use if "no" or undefined.

if(Backtrace_INCLUDE_DIR)
  # Already in cache, be silent
  set(Backtrace_FIND_QUIETLY TRUE)
endif(Backtrace_INCLUDE_DIR)

find_path(Backtrace_INCLUDE_DIR NAMES backtrace.h)

find_library(Backtrace_LIBRARY NAMES libbacktrace.a libbacktrace.so)

# handle the QUIETLY and REQUIRED arguments and set Backtrace_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Backtrace DEFAULT_MSG Backtrace_LIBRARY Backtrace_INCLUDE_DIR)

if(Backtrace_FOUND AND NOT TARGET Backtrace::Backtrace)
  add_library(Backtrace::Backtrace UNKNOWN IMPORTED)
  set_target_properties(
    Backtrace::Backtrace
    PROPERTIES IMPORTED_LOCATION "${Backtrace_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${Backtrace_INCLUDE_DIR}"
    )
endif()

if(Backtrace_FOUND)
  baci_add_dependency(baci_all_enabled_external_dependencies Backtrace::Backtrace)
  baci_add_dependency(baci_all_enabled_external_dependencies ${CMAKE_DL_LIBS})
  message(STATUS "Backtrace include directory: ${Backtrace_INCLUDE_DIR}")
  message(STATUS "Backtrace library directory: ${Backtrace_LIBRARY}")
endif(Backtrace_FOUND)
