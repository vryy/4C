/*-----------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometric search.

\level 3

*/
/*-----------------------------------------------------------*/

#include "../drt_io/io_pstream.H"

#include "geometric_search_utils.H"

namespace GEOMETRICSEARCH::UTILS
{
  void PrintGeometricSearchDetails(const Epetra_Comm& comm, const UTILS::GeometricSearchInfo info)
  {
    const int numproc = comm.NumProc();
    const int myrank = comm.MyPID();

    std::vector<int> primitive_size(numproc, 0), my_primitive_size(numproc, 0);
    std::vector<int> predicate_size(numproc, 0), my_predicate_size(numproc, 0);
    std::vector<int> coupling_pair_size(numproc, 0), my_coupling_pair_size(numproc, 0);

    my_primitive_size[myrank] = info.primitive_size;
    my_predicate_size[myrank] = info.predicate_size;
    my_coupling_pair_size[myrank] = info.coupling_pair_size;

    comm.SumAll(&my_primitive_size[0], &primitive_size[0], numproc);
    comm.SumAll(&my_predicate_size[0], &predicate_size[0], numproc);
    comm.SumAll(&my_coupling_pair_size[0], &coupling_pair_size[0], numproc);

    if (myrank == 0)
    {
      IO::cout(IO::verbose) << "   Collision search:" << IO::endl;
      IO::cout(IO::verbose) << "   +-----+------------+------------+--------------+" << IO::endl;
      IO::cout(IO::verbose) << "   | PID | predicates | primitives |  found pairs |" << IO::endl;
      IO::cout(IO::verbose) << "   +-----+------------+------------+--------------+" << IO::endl;

      for (int npid = 0; npid < numproc; ++npid)
      {
        IO::cout(IO::verbose) << "   | " << std::setw(3) << npid << " | " << std::setw(10)
                              << primitive_size[npid] << " | " << std::setw(10)
                              << predicate_size[npid] << " | " << std::setw(12)
                              << coupling_pair_size[npid] << " | " << IO::endl;
        IO::cout(IO::verbose) << "   +-----+------------+------------+--------------+" << IO::endl;
      }
      IO::cout(IO::verbose) << IO::endl;
    }
  }
}  // namespace GEOMETRICSEARCH::UTILS