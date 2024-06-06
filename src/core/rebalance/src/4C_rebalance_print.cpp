/*-----------------------------------------------------------*/
/*! \file

\brief Utility functions for the partitioning/rebalancing

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_rebalance_print.hpp"

#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

void Core::Rebalance::UTILS::print_parallel_distribution(const Discret::Discretization& dis)
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
      Core::IO::cout(Core::IO::verbose) << "\n   discretization: " << dis.Name() << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose)
          << "   +-----+---------------+--------------+-----------------+----------------+"
          << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose)
          << "   | PID |  n_rownodes   | n_ghostnodes |  n_rowelements  |   n_ghostele   |"
          << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose)
          << "   +-----+---------------+--------------+-----------------+----------------+"
          << Core::IO::endl;

      for (int npid = 0; npid < numproc; ++npid)
      {
        Core::IO::cout(Core::IO::verbose)
            << "   | " << std::setw(3) << npid << " | " << std::setw(13) << n_nodes[npid] << " | "
            << std::setw(12) << n_ghostnodes[npid] << " | " << std::setw(15) << n_elements[npid]
            << " | " << std::setw(14) << n_ghostele[npid] << " | " << Core::IO::endl;
        Core::IO::cout(Core::IO::verbose)
            << "   +-----+---------------+--------------+-----------------+----------------+"
            << Core::IO::endl;
      }
      Core::IO::cout(Core::IO::verbose) << Core::IO::endl;
    }
  }
}
FOUR_C_NAMESPACE_CLOSE
