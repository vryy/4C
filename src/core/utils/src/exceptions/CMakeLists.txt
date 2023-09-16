baci_auto_define_module()

baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} config)

find_package(Backtrace)

target_link_libraries(${AUTO_DEFINED_MODULE_NAME} PUBLIC ${CMAKE_DL_LIBS})
if(${Backtrace_FOUND})
  target_compile_definitions(${AUTO_DEFINED_MODULE_NAME} PRIVATE BOOST_STACKTRACE_USE_BACKTRACE)
  target_link_libraries(${AUTO_DEFINED_MODULE_NAME} PUBLIC Backtrace::Backtrace)
else()
  message(STATUS "libbacktrace is not found. You will not have line numbers in the stacktrace.")
  target_compile_definitions(${AUTO_DEFINED_MODULE_NAME} PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
endif()
