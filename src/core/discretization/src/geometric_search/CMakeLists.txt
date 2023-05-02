set(OBJS_GEOMETRIC_SEARCH
    # cmake-format: sortable
    discretization_geometric_search.cpp
    discretization_geometric_search_params.cpp
    discretization_geometric_search_utils.cpp
    )

set(LIBS_GEOMETRIC_SEARCH
    # cmake-format: sortable
    io
    lib
    )

baci_add_library_to_module(
  MODULE
  discretization
  NAME
  geometric_search
  SOURCES
  ${OBJS_GEOMETRIC_SEARCH}
  DEPENDENCIES
  ${LIBS_GEOMETRIC_SEARCH}
  )
