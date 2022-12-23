set(OBJS_GEOMETRY
    # cmake-format: sortable
    integrationcell.cpp
    intersection_interfacepoint.cpp
    intersection_math.cpp
    intersection_service.cpp
    position_array.cpp
    searchtree.cpp
    searchtree_geometry_service.cpp
    searchtree_nearestobject.cpp
    )

baci_add_library(geometry ${OBJS_GEOMETRY})

set(LIBS_GEOMETRY
    # cmake-format: sortable
    lib
    linalg
    )

target_link_libraries(geometry PUBLIC ${LIBS_GEOMETRY})
