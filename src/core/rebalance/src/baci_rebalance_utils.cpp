/*-----------------------------------------------------------*/
/*! \file

\brief Utility functions for the partitioning/rebalancing

\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_rebalance_utils.hpp"

#include "baci_io_pstream.hpp"
#include "baci_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

void CORE::REBALANCE::UTILS::PrintParallelDistribution(const DRT::Discretization& dis)
{
  const int numproc = dis.Comm().NumProc();
  const int myrank = dis.Comm().MyPID();

  if (numproc > 1)
  {
    std::vector<int> my_n_nodes(numproc, 0);
    std::vector<int> n_nodes(numproc, 0);
    std::vector<int> my_n_ghostnodes(numproc, 0);
    std::vector<int> n_ghostnodes(numproc, 0);
    std::vector<int> my_n_elements(numproc, 0);
    std::vector<int> n_elements(numproc, 0);
    std::vector<int> my_n_ghostele(numproc, 0);
    std::vector<int> n_ghostele(numproc, 0);

    my_n_nodes[myrank] = dis.NumMyRowNodes();
    my_n_ghostnodes[myrank] = dis.NumMyColNodes() - my_n_nodes[myrank];
    my_n_elements[myrank] = dis.NumMyRowElements();
    my_n_ghostele[myrank] = dis.NumMyColElements() - my_n_elements[myrank];

    dis.Comm().SumAll(&my_n_nodes[0], &n_nodes[0], numproc);
    dis.Comm().SumAll(&my_n_ghostnodes[0], &n_ghostnodes[0], numproc);
    dis.Comm().SumAll(&my_n_elements[0], &n_elements[0], numproc);
    dis.Comm().SumAll(&my_n_ghostele[0], &n_ghostele[0], numproc);

    if (myrank == 0)
    {
      IO::cout(IO::verbose) << "\n   Discretization: " << dis.Name() << IO::endl;
      IO::cout(IO::verbose)
          << "   +-----+---------------+--------------+-----------------+----------------+"
          << IO::endl;
      IO::cout(IO::verbose)
          << "   | PID |  n_rownodes   | n_ghostnodes |  n_rowelements  |   n_ghostele   |"
          << IO::endl;
      IO::cout(IO::verbose)
          << "   +-----+---------------+--------------+-----------------+----------------+"
          << IO::endl;

      for (int npid = 0; npid < numproc; ++npid)
      {
        IO::cout(IO::verbose) << "   | " << std::setw(3) << npid << " | " << std::setw(13)
                              << n_nodes[npid] << " | " << std::setw(12) << n_ghostnodes[npid]
                              << " | " << std::setw(15) << n_elements[npid] << " | "
                              << std::setw(14) << n_ghostele[npid] << " | " << IO::endl;
        IO::cout(IO::verbose)
            << "   +-----+---------------+--------------+-----------------+----------------+"
            << IO::endl;
      }
      IO::cout(IO::verbose) << IO::endl;
    }
  }
}
FOUR_C_NAMESPACE_CLOSE
