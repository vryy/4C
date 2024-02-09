set(OBJS_CREATERTD_MAIN
    # cmake-format: sortable
    baci_create_rtdfiles_base.cpp
    )
set(OBJS_CREATERTD_LIB baci_create_rtdfiles_wrapper.cpp baci_create_rtdfiles_utils.cpp)
baci_add_library(
  create_rtdfiles
  SOURCES
  ${OBJS_CREATERTD_LIB}
  HEADERS
  baci_create_rtdfiles_utils.hpp
  baci_create_rtdfiles_wrapper.hpp
  )

baci_add_dependency(create_rtdfiles config lib)
# the create_rtdfiles executable
add_executable(create_rtd EXCLUDE_FROM_ALL ${OBJS_CREATERTD_MAIN})
target_link_libraries(create_rtd PRIVATE create_rtdfiles_objs)
set_target_properties(create_rtd PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_target_requires_cyclic_legacy_linking(create_rtd)
add_custom_target(rtd DEPENDS create_rtd)
