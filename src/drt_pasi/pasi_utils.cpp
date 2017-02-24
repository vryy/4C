/*!----------------------------------------------------------------------
\file pasi_utils.cpp

\brief utility methods for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "pasi_utils.H"

/*----------------------------------------------------------------------*
 | modification of time parameter list                   sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::UTILS::ChangeTimeParameter(
    const Epetra_Comm& comm,                    //! local epetra communicator
    const Teuchos::ParameterList& pasi_params,  //! input parameters for particle structure interaction
    Teuchos::ParameterList& particle_params,    //! input parameters for particle dynamics
    Teuchos::ParameterList& struct_params       //! input parameters for structural dynamics
    )
{
  // set time parameters from particle structure interaction to subproblems

  // the default time step size
  particle_params.set<double>("TIMESTEP",pasi_params.get<double>("TIMESTEP"));
  struct_params.set<double>("TIMESTEP",pasi_params.get<double>("TIMESTEP"));

  // maximum number of timesteps
  particle_params.set<int>("NUMSTEP",pasi_params.get<int>("NUMSTEP"));
  struct_params.set<int>("NUMSTEP",pasi_params.get<int>("NUMSTEP"));

  // maximum simulation time
  particle_params.set<double>("MAXTIME",pasi_params.get<double>("MAXTIME"));
  struct_params.set<double>("MAXTIME",pasi_params.get<double>("MAXTIME"));

  // solution output
  particle_params.set<int>("RESULTSEVRY",pasi_params.get<int>("RESULTSEVRY"));
  struct_params.set<int>("RESULTSEVRY",pasi_params.get<int>("RESULTSEVRY"));

  // restart
  particle_params.set<int>("RESTARTEVRY",pasi_params.get<int>("RESTARTEVRY"));
  struct_params.set<int>("RESTARTEVRY",pasi_params.get<int>("RESTARTEVRY"));

  if (comm.MyPID() == 0)
  {
    std::cout << "================= Overview of chosen time stepping: ==================" << std::endl;
    std::cout << std::setw(20) << "Timestep:" << std::setw(15) << std::scientific  << std::setprecision(4) << pasi_params.get<double>("TIMESTEP") << std::endl;
    std::cout << std::setw(20) << "Numstep:" << std::setw(15) << std::scientific  << std::setprecision(4) << pasi_params.get<int>("NUMSTEP") << std::endl;
    std::cout << std::setw(20) << "Maxtime:" << std::setw(15) << std::scientific  << std::setprecision(4) << pasi_params.get<double>("MAXTIME") << std::endl;
    std::cout << std::setw(20) << "Result every step:" << std::setw(15) << pasi_params.get<int>("RESULTSEVRY") << std::endl;
    std::cout << std::setw(20) << "Restart every step:" << std::setw(15) << pasi_params.get<int>("RESTARTEVRY") << std::endl;
    std::cout << "======= currently equal for both structure and particle field ========" << std::endl;
  }

  return;
} // PASI::UTILS::ChangeTimeParameter()

/*----------------------------------------------------------------------*
 | print particle structure interaction logo             sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::UTILS::Logo()
{
  //FUCHS todo a Logo
  std::cout << "============================ Welcome to ==============================" << std::endl;
  std::cout << "  ___          _   _    _       ___ _               _                 " << std::endl;
  std::cout << " | _ \\__ _ _ _| |_(_)__| |___  / __| |_ _ _ _  _ __| |_ _  _ _ _ ___  " << std::endl;
  std::cout << " |  _/ _` | '_|  _| / _| / -_) \\__ \\  _| '_| || / _|  _| || | '_/ -_) " << std::endl;
  std::cout << " |_| \\__,_|_|  \\__|_\\__|_\\___| |___/\\__|_|  \\_,_\\__|\\__|\\_,_|_| \\___| " << std::endl;
  std::cout << "              ___     _                   _   _                       " << std::endl;
  std::cout << "             |_ _|_ _| |_ ___ _ _ __ _ __| |_(_)___ _ _               " << std::endl;
  std::cout << "              | || ' \\  _/ -_) '_/ _` / _|  _| / _ \\ ' \\              " << std::endl;
  std::cout << "             |___|_||_\\__\\___|_| \\__,_\\__|\\__|_\\___/_||_|             " << std::endl;
  std::cout << "" << std::endl;
  std::cout << "======================================================================" << std::endl;

  return;
} // PASI::UTILS::Logo()
