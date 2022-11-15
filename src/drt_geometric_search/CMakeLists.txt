set(OBJS_DRT_GEOMETRIC_SEARCH
    geometric_search.cpp geometric_search_params.cpp geometric_search_utils.cpp
    )

baci_add_library(drt_lib_geometric_search ${OBJS_DRT_GEOMETRIC_SEARCH})
target_link_libraries(drt_lib_geometric_search PUBLIC drt_lib)
