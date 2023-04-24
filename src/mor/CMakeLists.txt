set(OBJS_MOR
    # cmake-format: sortable
    mor_pod.cpp
    )

baci_add_library(mor ${OBJS_MOR})

set(LIBS_MOR
    # cmake-format: sortable
    linalg
    )

baci_add_dependency(mor ${LIBS_MOR})
