/*--------------------------------------------------------------------------*/
/*!

\brief utility class for  elastohydrodynamic lubrication (lubrication structure interaction)

\level 3

\maintainer Mostafa Faraji

*/
/*--------------------------------------------------------------------------*/

#include "ehl_utils.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 09/2014 */
/* Function for checking that the different time steps are a
 multiplicative of each other                                           */

int EHL::Utils::CheckTimeStepping(double dt1, double dt2)
{
  double workdt1 = std::min(dt1, dt2);
  double workdt2 = std::max(dt1, dt2);
  double t1 = 0.0;
  int i = 0;

  while (true)
  {
    i++;
    t1 = i * workdt1;

    if (std::abs(t1 - workdt2) < 10E-10)
      break;

    else if (t1 > workdt2)
      dserror("Chosen time steps %f and %f are not a multiplicative of each other", dt1, dt2);
  }
  return i;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 10/2014 */
// Modification of time parameter list for problem with different time step size

void EHL::Utils::ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ehlparams,
    Teuchos::ParameterList& lubricationdyn, Teuchos::ParameterList& sdyn)
{
  bool difftimestep = DRT::INPUT::IntegralValue<int>(ehlparams, "DIFFTIMESTEPSIZE");

  if (difftimestep)  // Create subproblems with different time steps
  {
    // Check correct choice of time stepping for single fields
    double lubricationstep = lubricationdyn.get<double>("TIMESTEP");
    double solidstep = sdyn.get<double>("TIMESTEP");

    CheckTimeStepping(lubricationstep, solidstep);

    // modify global time step size
    ehlparams.set<double>("TIMESTEP", std::min(lubricationstep, solidstep));
  }
  else
  {
    // -------------------------------------------------------------------
    // overrule certain parameters for coupled problems
    // -------------------------------------------------------------------
    // the default time step size
    lubricationdyn.set<double>("TIMESTEP", ehlparams.get<double>("TIMESTEP"));
    sdyn.set<double>("TIMESTEP", ehlparams.get<double>("TIMESTEP"));
    // maximum simulation time
    lubricationdyn.set<double>("MAXTIME", ehlparams.get<double>("MAXTIME"));
    sdyn.set<double>("MAXTIME", ehlparams.get<double>("MAXTIME"));
    // maximum number of timesteps
    lubricationdyn.set<int>("NUMSTEP", ehlparams.get<int>("NUMSTEP"));
    sdyn.set<int>("NUMSTEP", ehlparams.get<int>("NUMSTEP"));
  }

  // Check correct input of restart. Code relies that both time value RESTARTEVRYTIME and
  // RESULTSEVRYTIME are given if restart from time is applied
  double restarttime = ehlparams.get<double>("RESTARTEVRYTIME");
  double updatetime = ehlparams.get<double>("RESULTSEVRYTIME");
  if ((updatetime > 0.0) or (restarttime > 0.0))
    if (!(updatetime > 0.0) and !(restarttime > 0.0))
      dserror(
          "If time controlled output and restart is desired, both parameters RESTARTEVRYTIME and "
          "RESULTSEVRYTIME has to be set");

  // set restart params
  int lubricationrestart;
  int structurerestart;

  if (restarttime > 0.0)
  {
    lubricationrestart = CheckTimeStepping(lubricationdyn.get<double>("TIMESTEP"), restarttime);
    structurerestart = CheckTimeStepping(sdyn.get<double>("TIMESTEP"), restarttime);
  }
  else
  {
    int restart = ehlparams.get<int>("RESTARTEVRY");
    lubricationrestart = restart;
    structurerestart = restart;
  }

  // set output params
  int lubricationupres;
  int structureupres;

  if (updatetime > 0.0)
  {
    lubricationupres = CheckTimeStepping(lubricationdyn.get<double>("TIMESTEP"), updatetime);
    structureupres = CheckTimeStepping(sdyn.get<double>("TIMESTEP"), updatetime);
  }
  else
  {
    int update = ehlparams.get<int>("RESULTSEVRY");
    lubricationupres = update;
    structureupres = update;
  }

  // restart
  lubricationdyn.set<int>("RESTARTEVRY", lubricationrestart);
  sdyn.set<int>("RESTARTEVRY", structurerestart);
  // solution output
  lubricationdyn.set<int>("RESULTSEVRY", lubricationupres);
  sdyn.set<int>("RESULTSEVRY", structureupres);

  if (comm.MyPID() == 0)
  {
    std::cout << "====================== Overview of chosen time stepping: "
                 "==============================\n"
              << "\t Timestep lubrication:           " << lubricationdyn.get<double>("TIMESTEP")
              << "\n"
              << "\t Timestep structure:        " << sdyn.get<double>("TIMESTEP") << "\n"
              << "\t Result step lubrication:        " << lubricationdyn.get<int>("RESULTSEVRY")
              << "\n"
              << "\t Result step structure:     " << sdyn.get<int>("RESULTSEVRY") << "\n"
              << "\t Restart step lubrication:       " << lubricationdyn.get<int>("RESTARTEVRY")
              << "\n"
              << "\t Restart step structure:    " << sdyn.get<int>("RESTARTEVRY") << "\n"
              << "================================================================================="
                 "=======\n \n";
  }
}

/*----------------------------------------------------------------------*
 | print EHL-logo                                            Faraji 05/19 |
 *----------------------------------------------------------------------*/
void EHL::printlogo()
{
  // more at http://www.ascii-art.de
  std::cout
      << "                                    ,/                                              \n"
      << "                                    #%                                              \n"
      << "                                   ,&&                                              \n"
      << "                                  ,&%%(                                             \n"
      << "                                  %&&&&%                                            \n"
      << "                                 /&&%%%%#                                           \n"
      << "                                ,&&%&&&&%#                                          \n"
      << "                                %&&&&&&&&&.                                         \n"
      << "                               .#&&&&&&&%%*                                         \n"
      << "                                %&&%&&&&%&                                          \n"
      << "                     *%%(        *&&%%%&/           .&&&%&                          \n"
      << "            .#&%,    &&&&%,                 /*(     *%&&&&,    %&&&%                \n"
      << "           %%%&&%%%&&&&&&&,  /&&&          %&&&&/  /%%&&&&&&%#&&&&&&                \n"
      << "            #%&&&&%%&%&&&&&#&%&&&          &&&&&%&&&&&&&&&&&&&&&&/&&                \n"
      << "           .#&&&%,   *#&&&&&&%&&/          %&&&&&&&&&&&&&&&&&&&&&*%&                \n"
      << "       &%&%&&&%%.         /&&&&%%*   *&&  ,%&&&&&&%&&&%###%&&%%&&&&&&%%&&&%&        \n"
      << "      *%&&%&&%%            .&&&&.   /&&&&%&&&&%&&&(           .%&&&&&&&&&&&/        \n"
      << "      (&&%&&&&*             %%&%%,  %&&&&&&&&&&&/                #&&&&&&&/          \n"
      << "         .%&&&*             %%&&&&&&. *%%&&&&&&                   (&&&&&&&          \n"
      << "          /%&%&            *&&&&&&&&   %&&&&&&.                    /%&&&&&&&&       \n"
      << "        ,&&&&&&&(        .%&&&&(,**,.*#&&&&&%%                      &&&&&&&&&&&.    \n"
      << "       .%%&%&&&&&%&%#(#%&&&&&&    #%&&&&&&&&%%                      &%&&&&&%&&%,    \n"
      << "         #&&( ,&&&&&&%%&&&&&&%(   /&&&&&&&&&%&                     .&&&&&%/         \n"
      << "                (%&&&%*  ,&&&.       /%&&&&%(                   ,%&&&&&&,           \n"
      << "                %&%%&,     *.           %&&&&&&(                 ,&&&&&&&&%&(.      \n"
      << "                                      .%&&&&&&&&&(             ,%%&&&&&&&&&%&/      \n"
      << "                                     (%&&&&&&&&&&&&&%/,    *(%&&&&&&&%&&&%&&,       \n"
      << "                                      (&&&%(%&&&&&&&&&%&&&&&&&&&&&&&                \n"
      << "                                              %&&&&&&&&&&&&&&&&&&&&&%&              \n"
      << "                                             .&&&&%&&%%&&&&&&&&(&&%&&&%.            \n"
      << "                                             #%&%%%    %%&&&&     %,                \n"
      << "                                               ,(*     *&&&&%                       \n"
      << "                                                        ...                         \n"
      << "                                                                                    \n"
      << "\n"
      << std::endl;

}  // printlogo()

/*----------------------------------------------------------------------*/
