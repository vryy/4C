/*---------------------------------------------------------------------*/
/*! \file

\brief Utility to write restart information based on a wall time interval

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_IO_WALLTIME_BASED_RESTART_HPP
#define FOUR_C_IO_WALLTIME_BASED_RESTART_HPP

#include "4C_config.hpp"

#include <Epetra_Comm.h>
#include <signal.h>
#include <stdio.h>

#include <random>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SparseMatrix;
}

namespace CORE::IO
{
  /*!
  \brief handles restart after a certain walltime interval, step interval or on a user signal

  \author hammerl
  */
  class RestartManager
  {
   public:
    RestartManager();

    virtual ~RestartManager() = default;

    /// setup of restart manager
    void SetupRestartManager(const double restartinterval, const int restartevry);

    /// return whether it is time for a restart
    /// \param step [in] : current time step for multi-field syncronisation
    /// \param comm [in] : get access to involved procs
    bool Restart(const int step, const Epetra_Comm& comm);

    /// the signal handler that gets passed to the kernel and listens for SIGUSR1 and SIGUSR2
    static void restart_signal_handler(
        int signal_number, siginfo_t* signal_information, void* ignored);

   private:
    /// @name wall time parameters
    //@{

    /// start time of simulation
    double startwalltime_;

    /// after this wall time interval a restart is enforced
    double restartevrytime_;

    /// check to enforce restart only once during time interval
    int restartcounter_;

    //@}

    /// store the step which was allowed to write restart
    int lastacceptedstep_;

    /// member to detect time step increment
    int lasttestedstep_;

    /// after this number of steps a restart is enforced
    int restartevrystep_;

    /// signal which was caught by the signal handler
    volatile static int signal_;
  };
}  // namespace CORE::IO


FOUR_C_NAMESPACE_CLOSE

#endif
