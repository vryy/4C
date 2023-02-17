set(OBJS_SURFSTRESS
    # cmake-format: sortable
    surfstress_manager.cpp
    )

baci_add_library(lib_surfstress ${OBJS_SURFSTRESS})

set(LIBS_LIB_STURFSTRESS
    # cmake-format: sortable
    lib
    )

baci_add_dependency(lib_surfstress ${LIBS_LIB_STURFSTRESS})
