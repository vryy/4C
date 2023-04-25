set(OBJS_MODULE_REGISTRY
    # cmake-format: sortable
    module_registry_parobjectregister.cpp
    )

baci_add_library(module_registry ${OBJS_MODULE_REGISTRY})

set(LIBS_MODULE_REGISTRY
    # cmake-format: sortable
    lib
    lib_constraint
    )

target_link_libraries(module_registry PUBLIC ${LIBS_MODULE_REGISTRY})
