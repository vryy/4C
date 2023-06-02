set(OBJS_CREATERTD_MAIN
    # cmake-format: sortable
    create_rtdfiles_base.cpp
    )
set(OBJS_CREATERTD_LIB create_rtdfiles_wrapper.cpp create_rtdfiles_utils.cpp)
baci_add_library(
  create_rtdfiles
  SOURCES
  ${OBJS_CREATERTD_LIB}
  HEADERS
  create_rtdfiles_utils.H
  create_rtdfiles_wrapper.H
  )

baci_add_dependency(create_rtdfiles config lib)
# the create_rtdfiles executable
add_executable(create_rtd EXCLUDE_FROM_ALL ${OBJS_CREATERTD_MAIN})
baci_add_dependency(create_rtd create_rtdfiles)
set_target_properties(create_rtd PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_target_requires_cyclic_legacy_linking(create_rtd)
add_custom_target(rtd DEPENDS create_rtd)
