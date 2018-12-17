add_library(geometry_pair_unittest_stubs
DatFileReader_stub.cpp
)
target_compile_options( geometry_pair_unittest_stubs PRIVATE -g -O0 )
target_include_directories( geometry_pair_unittest_stubs PRIVATE ${PROJECT_SOURCE_DIR} )

add_library(geometry_pair_dserror_stub dserror_stub.cpp)
target_compile_options( geometry_pair_dserror_stub PRIVATE -g -O0 )
target_include_directories( geometry_pair_dserror_stub PRIVATE ${PROJECT_SOURCE_DIR} )
