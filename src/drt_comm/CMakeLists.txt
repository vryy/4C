set(OBJS_DRT_COMM comm_utils.cpp)

add_library(drt_comm ${OBJS_DRT_COMM})
baci_link_default_libraries(drt_comm)
target_include_directories(drt_comm PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(drt_comm PUBLIC drt_lib)
