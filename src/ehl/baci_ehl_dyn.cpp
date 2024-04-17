/*--------------------------------------------------------------------------*/
/*! \file

\brief control routine for elastohydrodynamic lubrication (lubrication structure interaction)

\level 3

*/
/*--------------------------------------------------------------------------*/

#include "baci_ehl_dyn.hpp"

#include "baci_ehl_monolithic.hpp"
#include "baci_ehl_partitioned.hpp"
#include "baci_ehl_utils.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | BACI Logo for  EHL problems                             Faraji 05/19 |
 *----------------------------------------------------------------------*/
void printehllogo()
{
  std::cout << "---------------------------------------------------------------------------------"
            << std::endl;
  std::cout << "---------------------------------------------------------------------------------"
            << std::endl;
  std::cout << "-----------  Welcome to the Elasto-Hydrodynamic Lubrication problem  ------------"
            << std::endl;
  std::cout << "---------------------------------------------------------------------------------"
            << std::endl;
  std::cout << "---------------------------------------------------------------------------------"
            << std::endl;
  return;
}

void printehlmixlogo()
{
  std::cout << "---------------------------------------------------------------------------------"
            << std::endl;
  std::cout << "-----------------        Welcome to the problem type EHL        -----------------"
            << std::endl;
  std::cout << "-----------------               Mixed Lubrication               -----------------"
            << std::endl;
  std::cout << "-----------------           Averaged Reynolds Equation          -----------------"
            << std::endl;
  std::cout << "---------------------------------------------------------------------------------"
            << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Main control routine for EHL problems                    wirtz 12/15 |
 *----------------------------------------------------------------------*/
void ehl_dyn()
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // 1.- Initialization
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // 2.- Parameter reading
  Teuchos::ParameterList& ehlparams =
      const_cast<Teuchos::ParameterList&>(problem->ElastoHydroDynamicParams());
  // access lubrication params list
  Teuchos::ParameterList& lubricationdyn =
      const_cast<Teuchos::ParameterList&>(problem->LubricationDynamicParams());
  // do we want to use Modified Reynolds Equation?
  bool modifiedreynolds =
      (CORE::UTILS::IntegralValue<int>(lubricationdyn, "MODIFIED_REYNOLDS_EQU"));

  // print problem specific logo
  if (!problem->GetDis("structure")->Comm().MyPID())
  {
    if (!modifiedreynolds)
      printehllogo();
    else
      printehlmixlogo();
  }

  if (!problem->GetDis("structure")->Comm().MyPID()) EHL::printlogo();

  // access structural dynamic params list which will be possibly modified while creating the time
  // integrator
  Teuchos::ParameterList& sdyn =
      const_cast<Teuchos::ParameterList&>(GLOBAL::Problem::Instance()->StructuralDynamicParams());


  //  //Modification of time parameter list
  EHL::Utils::ChangeTimeParameter(comm, ehlparams, lubricationdyn, sdyn);

  const INPAR::EHL::SolutionSchemeOverFields coupling =
      CORE::UTILS::IntegralValue<INPAR::EHL::SolutionSchemeOverFields>(ehlparams, "COUPALGO");

  // 3.- Creation of Lubrication + Structure problem. (Discretization called inside)
  Teuchos::RCP<EHL::Base> ehl = Teuchos::null;

  // 3.1 choose algorithm depending on solution type
  switch (coupling)
  {
    case INPAR::EHL::ehl_IterStagg:
      ehl = Teuchos::rcp(
          new EHL::Partitioned(comm, ehlparams, lubricationdyn, sdyn, "structure", "lubrication"));
      break;
    case INPAR::EHL::ehl_Monolithic:
      ehl = Teuchos::rcp(
          new EHL::Monolithic(comm, ehlparams, lubricationdyn, sdyn, "structure", "lubrication"));
      break;
    default:
      dserror("unknown coupling algorithm for EHL!");
      break;
  }

  // 3.2- Read restart if needed. (Discretization called inside)
  const int restart = problem->Restart();
  if (restart) ehl->ReadRestart(restart);

  // 4.- Run of the actual problem.

  // 4.1.- Some setup needed for the elastohydrodynamic lubrication problem.
  ehl->SetupSystem();

  // 4.2.- Solve the whole problem
  ehl->Timeloop();

  // 4.3.- Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  ehl->TestResults(comm);
}

FOUR_C_NAMESPACE_CLOSE
