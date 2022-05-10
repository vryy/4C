#include <cxxtest/RealDescriptions.h>
#include <cxxtest/TestMain.h>
#include <cxxtest/ErrorPrinter.h>
#include <mpi.h>

int main( int argc, char *argv[] )
{
  int status;
  CxxTest::ErrorPrinter tmp;
  CxxTest::RealWorldDescription::_worldName = "cxxtest";
  MPI_Init(&argc, &argv);
  status = CxxTest::Main< CxxTest::ErrorPrinter >( tmp, argc, argv );
  MPI_Finalize();
  return status;
}

<CxxTest world>
