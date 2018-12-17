# This is the base of all unittests
find_package( CxxTest 4 )
if ( CXXTEST_FOUND )
  enable_testing()
  add_subdirectory ( Unittests )
  add_subdirectory ( tests/geometry_pair_tests )
else ()
  message(" WW Warning: Unittests are unavailable, as cxxtest is missing. Try: sudo yum install cxxtest")
endif ( CXXTEST_FOUND )