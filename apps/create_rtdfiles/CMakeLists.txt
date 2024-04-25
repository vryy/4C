set(OBJS_CREATERTD_MAIN
    # cmake-format: sortable
    4C_create_rtdfiles_base.cpp
    )
set(OBJS_CREATERTD_LIB 4C_create_rtdfiles_wrapper.cpp 4C_create_rtdfiles_utils.cpp)
baci_add_library(
  create_rtdfiles
  SOURCES
  ${OBJS_CREATERTD_LIB}
  HEADERS
  4C_create_rtdfiles_utils.hpp
  4C_create_rtdfiles_wrapper.hpp
  )

baci_add_dependency(create_rtdfiles config lib)
# the create_rtdfiles executable
add_executable(create_rtd EXCLUDE_FROM_ALL ${OBJS_CREATERTD_MAIN})
target_link_libraries(create_rtd PRIVATE create_rtdfiles_objs)
set_target_properties(create_rtd PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_set_up_executable(create_rtd)
add_custom_target(rtd DEPENDS create_rtd)
