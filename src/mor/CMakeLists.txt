set(OBJS_MOR
    # cmake-format: sortable
    mor_pod.cpp
    )

baci_add_library(mor ${OBJS_MOR})

set(LIBS_MOR
    # cmake-format: sortable
    linalg
    )

target_link_libraries(mor PUBLIC ${LIBS_MOR})
