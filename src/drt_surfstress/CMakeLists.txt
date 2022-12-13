set(OBJS_DRT_SURFSTRESS
    # cmake-format: sortable
    drt_surfstress_manager.cpp
    )

baci_add_library(drt_lib_surfstress ${OBJS_DRT_SURFSTRESS})

set(LIBS_DRT_LIB_STURFSTRESS
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(drt_lib_surfstress PUBLIC ${LIBS_DRT_LIB_STURFSTRESS})
