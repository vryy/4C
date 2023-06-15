baci_add_library(
  post_common
  SOURCES
  post_common.cpp
  HEADERS
  post_common.H
  )
baci_add_dependency(post_common inpar module_registry)

set(OBJS_POST_PROCESSOR_MAIN
    # cmake-format: sortable
    post_common_filter_base.cpp
    post_common_post_processor.cpp
    post_common_single_field_writers.cpp
    post_common_structure_stress.cpp
    post_common_thermo_heatflux.cpp
    post_common_writer_base.cpp
    )

add_executable(post_processor EXCLUDE_FROM_ALL ${OBJS_POST_PROCESSOR_MAIN})
set_target_properties(post_processor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_add_dependency(post_processor post_ensight_lib post_vtk_lib post_common)
baci_target_requires_cyclic_legacy_linking(post_processor)

# various post-processing tools
include(scripts/baci_create_post_scripts.cmake)
