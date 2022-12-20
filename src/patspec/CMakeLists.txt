set(OBJS_PATSPEC
    # cmake-format: sortable
    patspec.cpp
    )

baci_add_library(patspec ${OBJS_PATSPEC})

set(LIBS_PATSPEC
    # cmake-format: sortable
    mat
    )

target_link_libraries(patspec PUBLIC ${LIBS_PATSPEC})
