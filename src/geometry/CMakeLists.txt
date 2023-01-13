set(OBJS_GEOMETRY
    # cmake-format: sortable
    geometry_integrationcell.cpp
    geometry_intersection_interfacepoint.cpp
    geometry_intersection_math.cpp
    geometry_intersection_service.cpp
    geometry_position_array.cpp
    geometry_searchtree.cpp
    geometry_searchtree_nearestobject.cpp
    geometry_searchtree_service.cpp
    )

baci_add_library(geometry ${OBJS_GEOMETRY})

set(LIBS_GEOMETRY
    # cmake-format: sortable
    lib
    linalg
    )

target_link_libraries(geometry PUBLIC ${LIBS_GEOMETRY})
