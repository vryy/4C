find_package(Backtrace REQUIRED)

if(Backtrace_FOUND)
  message(STATUS "Backtrace include directory: ${Backtrace_INCLUDE_DIR}")
  message(STATUS "Backtrace library directory: ${Backtrace_LIBRARY}")

  target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE Backtrace::Backtrace)
endif()
