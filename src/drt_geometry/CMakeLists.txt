set(OBJS_DRT_GEOMETRY
    integrationcell.cpp
    intersection_interfacepoint.cpp
    intersection_math.cpp
    intersection_service.cpp
    position_array.cpp
    searchtree.cpp
    searchtree_geometry_service.cpp
    searchtree_nearestobject.cpp
    )

baci_add_library(drt_geometry ${OBJS_DRT_GEOMETRY})
target_link_libraries(drt_geometry PUBLIC linalg drt_lib)
