add_library(unittest_stubs
DatFileReader_stub.cpp
)
target_compile_options( unittest_stubs PRIVATE -g -O0 )
target_include_directories( unittest_stubs PRIVATE ${PROJECT_SOURCE_DIR} )

add_library(dserror_stub dserror_stub.cpp)
target_compile_options( dserror_stub PRIVATE -g -O0 )
target_include_directories( dserror_stub PRIVATE ${PROJECT_SOURCE_DIR} )
