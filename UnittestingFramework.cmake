# This is the base of all unittests
find_package( CxxTest 4 )
if ( CxxTest_FOUND )
  enable_testing()
  add_subdirectory ( Unittests )
else ()
  message(" WW Warning: Unittests are unavailable, as cxxtest is missing. Try: sudo dnf install cxxtest")
endif ( CxxTest_FOUND )

