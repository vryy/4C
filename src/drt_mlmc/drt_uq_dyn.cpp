/*----------------------------------------------------------------------*/
/*! \file
\brief Uncertainty Quantification

\maintainer Jonas Nitzler

\level 3
*/
/*----------------------------------------------------------------------*/

#ifdef HAVE_FFTW


/*----------------------------------------------------------------------*
 * Main control routine for uq here we decide what problem we want to solve
 *
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "drt_uq_dyn.H"
#include "../drt_lib/drt_globalproblem.H"
#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io.H"
#include "../drt_red_airways/airwayimplicitintegration.H"

#include "../drt_red_airways/red_airway_resulttest.H"
#include "../drt_red_airways/red_airways_dyn_drt.H"
#include "../drt_red_airways/airwayimplicitintegration.H"
#include "../drt_red_airways/redairway_tissue.H"
#include "../drt_adapter/ad_str_redairway.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_mlmc/drt_uq_redairways.H"
#include "../drt_mlmc/mlmc.H"
#include <ctime>
#include <cstdlib>
#include <iostream>


void dyn_uq()
{
  // get list for multi level monte carlo
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  INPAR::MLMC::FWDProblem fwdprb =
      DRT::INPUT::IntegralValue<INPAR::MLMC::FWDProblem>(mlmcp, "FWDPROBLEM");

  // quick check whether to read in structure or airways
  if (fwdprb == INPAR::MLMC::structure)
  {
    IO::cout << "setup " << IO::endl;
    structure_uq();
  }
  else if (fwdprb == INPAR::MLMC::red_airways)
  {
    redairways_uq();
  }
  else
    dserror("Unknown forward problem type fix your input file");
}


void redairways_uq()
{
  if (DRT::Problem::Instance()->DoesExistDis("red_airway") == false)
  {
    dserror("Red Airways do not exist");
    // return Teuchos::null;
  }
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("red_airway");

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled())
  {
    actdis->FillComplete();
  }

  // -------------------------------------------------------------------
  // If discretization is empty, then return empty time integration
  // -------------------------------------------------------------------
  if (actdis->NumGlobalElements() < 1)
  {
    dserror("Red Airways Dis is empty");
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  // Teuchos::RCP<IO::DiscretizationWriter>  output = actdis->Writer();
  // output->WriteMesh(0,0.0);
  // get input lists
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();


  INPAR::MLMC::UQStrategy uqstrat =
      DRT::INPUT::IntegralValue<INPAR::MLMC::UQStrategy>(mlmcp, "UQSTRATEGY");
  UQ::UQ_REDAIRWAYS mc(actdis);

  switch (uqstrat)
  {
    case INPAR::MLMC::mc_plain:
    {
      mc.Integrate();
    }
    break;
    case INPAR::MLMC::mc_paracont:
    {
      dserror("Unknown UQ Strategy for redairways problem fix your input file ");
    }
    break;
    case INPAR::MLMC::mc_scaledthick:
    {
      dserror("Unknown UQ Strategy for redairways problem fix your input file ");
    }
    break;
    default:
      dserror("Unknown UQ Strategy for redairways problem fix your input file ");
      break;
  }

}  // end of redairways_uq()


/*======================================================================*/
/* Monte Carlo analysis of structures */
void structure_uq()
{
  // get input lists
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("structure");

  // set degrees of freedom in the discretization
  if (not actdis->Filled() || not actdis->HaveDofs()) actdis->FillComplete();

  // input parameters for structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // show default parameters
  if (actdis->Comm().MyPID() == 0) DRT::INPUT::PrintDefaultParameters(IO::cout, sdyn);

  // create a solver
  // get the solver number used for structural solver
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL "
        "DYNAMIC to a valid number!");

  INPAR::MLMC::UQStrategy uqstrat =
      DRT::INPUT::IntegralValue<INPAR::MLMC::UQStrategy>(mlmcp, "UQSTRATEGY");
  UQ::MLMC mc(actdis);

  switch (uqstrat)
  {
    case INPAR::MLMC::mc_plain:
    {
      mc.Integrate();
    }
    break;
    case INPAR::MLMC::mc_paracont:
    {
      mc.IntegrateNoReset();
    }
    break;
    case INPAR::MLMC::mc_scaledthick:
    {
      mc.IntegrateScaleByThickness();
    }
    break;
    default:
      dserror("Unknown UQ Strategy for Structure problem fix your input file ");
      break;
  }


}  // end of structure_uq()


#endif  // FFTW
