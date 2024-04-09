/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_lib_utils.hpp"

#include "baci_discretization_fem_general_element_center.hpp"
#include "baci_global_data.hpp"
#include "baci_io_control.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_FEVector.h>


BACI_NAMESPACE_OPEN


DRT::UTILS::RestartManager::RestartManager()
    : startwalltime_(GLOBAL::Problem::Walltime()),
      restartevrytime_(-1.0),
      restartcounter_(0),
      lastacceptedstep_(-1),
      lasttestedstep_(-1),
      restartevrystep_(-1)
{
  // setup signal handler
  signal_ = -1;
  struct sigaction the_action;
  the_action.sa_sigaction = restart_signal_handler;
  sigemptyset(&the_action.sa_mask);
  the_action.sa_flags = SA_SIGINFO;

  if (sigaction(SIGUSR1, &the_action, nullptr))
    dserror("signal handler for action SIGUSR1 could not be registered");
  if (sigaction(SIGUSR2, &the_action, nullptr))
    dserror("signal handler for action SIGUSR2 could not be registered");
}

/// set the time interval to enforce restart writing
void DRT::UTILS::RestartManager::SetupRestartManager(
    const double restartinterval, const int restartevry)
{
  restartevrytime_ = restartinterval;
  restartevrystep_ = restartevry;
}

/// return whether it is time for a restart after a certain walltime interval
bool DRT::UTILS::RestartManager::Restart(const int step, const Epetra_Comm& comm)
{
  // make sure that all after the first field write restart, too
  if (step == lastacceptedstep_) return true;

  // make sure that only the first field tests the time limit
  if (step > lasttestedstep_)
  {
    lasttestedstep_ = step;

    // compute elapsed walltime on proc 0 and let it decide for all other procs, too
    int restarttime = 0;
    if (comm.MyPID() == 0)
    {
      const double elapsedtime = GLOBAL::Problem::Walltime() - startwalltime_;
      const bool walltimerestart = (int)(elapsedtime / restartevrytime_) > restartcounter_;

      if (step > 0 and (((restartevrystep_ > 0) and (step % restartevrystep_ == 0)) or
                           walltimerestart or signal_ > 0))
      {
        lastacceptedstep_ = step;
        restarttime = 1;
        signal_ = -1;
        // only increment counter for walltime based restart functionality
        if (walltimerestart) ++restartcounter_;
      }
    }
    comm.Broadcast(&restarttime, 1, 0);
    return restarttime;
  }

  return false;
}

void DRT::UTILS::RestartManager::restart_signal_handler(
    int signal_number, siginfo_t* signal_information, void* ignored)
{
  signal_ = signal_information->si_signo;
  return;
}

volatile int DRT::UTILS::RestartManager::signal_;

BACI_NAMESPACE_CLOSE
