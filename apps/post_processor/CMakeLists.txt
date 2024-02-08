set(OBJS_POST_PROCESSOR_MAIN
    # cmake-format: sortable
    baci_post_processor.cpp
    baci_post_processor_single_field_writers.cpp
    baci_post_processor_structure_stress.cpp
    baci_post_processor_thermo_heatflux.cpp
    )

add_executable(post_processor EXCLUDE_FROM_ALL ${OBJS_POST_PROCESSOR_MAIN})
set_target_properties(post_processor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
baci_add_dependency(post_processor post)
baci_target_requires_cyclic_legacy_linking(post_processor)

# various post-processing tools
include(scripts/baci_create_post_scripts.cmake)
