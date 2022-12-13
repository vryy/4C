set(OBJS_DRT_GEOMETRIC_SEARCH
    # cmake-format: sortable
    geometric_search.cpp
    geometric_search_params.cpp
    geometric_search_utils.cpp
    )

baci_add_library(drt_lib_geometric_search ${OBJS_DRT_GEOMETRIC_SEARCH})

set(LIBS_DRT_LIB_GEOMETRIC_SEARCH
    # cmake-format: sortable
    drt_io
    drt_lib
    )

target_link_libraries(drt_lib_geometric_search PUBLIC ${LIBS_DRT_LIB_GEOMETRIC_SEARCH})
