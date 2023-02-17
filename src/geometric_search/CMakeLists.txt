set(OBJS_GEOMETRIC_SEARCH
    # cmake-format: sortable
    geometric_search.cpp
    geometric_search_params.cpp
    geometric_search_utils.cpp
    )

baci_add_library(lib_geometric_search ${OBJS_GEOMETRIC_SEARCH})

set(LIBS_LIB_GEOMETRIC_SEARCH
    # cmake-format: sortable
    io
    lib
    )

baci_add_dependency(lib_geometric_search ${LIBS_LIB_GEOMETRIC_SEARCH})
