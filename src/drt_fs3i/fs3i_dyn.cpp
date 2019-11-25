/*----------------------------------------------------------------------*/
/*! \file
\brief Main control routine for fluid-structure-scalar-scalar
       interaction (FS3I)

\level 1

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/


#include "fs3i_dyn.H"

#include "fs3i.H"
#include "fs3i_partitioned_1wc.H"
#include "fps3i_partitioned_1wc.H"
#include "fs3i_partitioned_2wc.H"
#include "biofilm_fsi.H"
#include "ac_fsi.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
// entry point for all kinds of FS3I
/*----------------------------------------------------------------------*/
void fs3i_dyn()
{
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("structure")->Comm();

  Teuchos::RCP<FS3I::FS3I_Base> fs3i;

  // what's the current problem type?
  ProblemType probtype = DRT::Problem::Instance()->GetProblemType();

  switch (probtype)
  {
    case prb_gas_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::PartFS3I_1WC(comm));
    }
    break;
    case prb_ac_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::ACFSI(comm));
    }
    break;
    case prb_thermo_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::PartFS3I_2WC(comm));
    }
    break;
    case prb_biofilm_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::BiofilmFSI(comm));
    }
    break;
    case prb_fps3i:
    {
      fs3i = Teuchos::rcp(new FS3I::PartFPS3I_1WC(comm));
    }
    break;
    default:
      dserror("solution of unknown problemtyp %d requested", probtype);
      break;
  }

  fs3i->Init();
  fs3i->Setup();

  // read the restart information, set vectors and variables ---
  // be careful, dofmaps might be changed here in a Redistribute call
  fs3i->ReadRestart();

  // if running FPS3I in parallel one needs to redistribute the interface after restarting
  fs3i->RedistributeInterface();

  // now do the coupling and create combined dofmaps
  fs3i->SetupSystem();

  fs3i->Timeloop();

  fs3i->TestResults(comm);

  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
}
