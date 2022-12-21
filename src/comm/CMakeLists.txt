set(OBJS_COMM
    # cmake-format: sortable
    comm_utils.cpp
    )

add_library(comm ${OBJS_COMM})
baci_link_default_libraries(comm)
target_include_directories(comm PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set(LIBS_COMM
    # cmake-format: sortable
    lib
    nurbs_discret
    rebalance
    )

target_link_libraries(comm PUBLIC ${LIBS_COMM})
