set(OBJS_DRT_MOR
    # cmake-format: sortable
    mor_pod.cpp
    )

baci_add_library(drt_mor ${OBJS_DRT_MOR})

set(LIBS_DRT_MOR
    # cmake-format: sortable
    linalg
    )

target_link_libraries(drt_mor PUBLIC ${LIBS_DRT_MOR})
