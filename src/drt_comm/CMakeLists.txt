set(OBJS_DRT_COMM
    # cmake-format: sortable
    comm_utils.cpp
    )

add_library(drt_comm ${OBJS_DRT_COMM})
baci_link_default_libraries(drt_comm)
target_include_directories(drt_comm PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set(LIBS_DRT_COMM
    # cmake-format: sortable
    drt_lib
    drt_nurbs_discret
    rebalance
    )

target_link_libraries(drt_comm PUBLIC ${LIBS_DRT_COMM})
