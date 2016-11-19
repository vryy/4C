# This is the base of all unittests
find_package( CxxTest 4 REQUIRED )
if ( CXXTEST_FOUND )
  enable_testing()
  add_subdirectory ( Unittests )
else ()
  message(" WW Warning: Unittests are unavailable, as cxxtest is missing. Try: sudo dnf install cxxtest")
endif ( CXXTEST_FOUND )

