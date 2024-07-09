/*---------------------------------------------------------------------*/
/*! \file

\brief Utility to write restart information based on a wall time interval

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_io_walltime_based_restart.hpp"

#include "4C_utils_exceptions.hpp"

#include <chrono>

FOUR_C_NAMESPACE_OPEN

namespace
{
  double walltime_in_seconds()
  {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::high_resolution_clock::now().time_since_epoch())
               .count() *
           1.0e-3;
  }
}  // namespace


Core::IO::RestartManager::RestartManager()
    : startwalltime_(walltime_in_seconds()),
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
    FOUR_C_THROW("signal handler for action SIGUSR1 could not be registered");
  if (sigaction(SIGUSR2, &the_action, nullptr))
    FOUR_C_THROW("signal handler for action SIGUSR2 could not be registered");
}

/// set the time interval to enforce restart writing
void Core::IO::RestartManager::setup_restart_manager(
    const double restartinterval, const int restartevry)
{
  restartevrytime_ = restartinterval;
  restartevrystep_ = restartevry;
}

/// return whether it is time for a restart after a certain walltime interval
bool Core::IO::RestartManager::restart(const int step, const Epetra_Comm& comm)
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
      const double elapsedtime = walltime_in_seconds() - startwalltime_;
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

void Core::IO::RestartManager::restart_signal_handler(
    int signal_number, siginfo_t* signal_information, void* ignored)
{
  signal_ = signal_information->si_signo;
}

volatile int Core::IO::RestartManager::signal_;

FOUR_C_NAMESPACE_CLOSE
