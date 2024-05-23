/*----------------------------------------------------------------------*/
/*! \file
\brief Main control routine for fluid-structure-scalar-scalar
       interaction (FS3I)

\level 1


*----------------------------------------------------------------------*/


#include "4C_fs3i_dyn.hpp"

#include "4C_comm_utils.hpp"
#include "4C_fs3i.hpp"
#include "4C_fs3i_ac_fsi.hpp"
#include "4C_fs3i_biofilm_fsi.hpp"
#include "4C_fs3i_fps3i_partitioned_1wc.hpp"
#include "4C_fs3i_partitioned_1wc.hpp"
#include "4C_fs3i_partitioned_2wc.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
// entry point for all kinds of FS3I
/*----------------------------------------------------------------------*/
void fs3i_dyn()
{
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("structure")->Comm();

  Teuchos::RCP<FS3I::FS3IBase> fs3i;

  // what's the current problem type?
  GLOBAL::ProblemType probtype = GLOBAL::Problem::Instance()->GetProblemType();

  switch (probtype)
  {
    case GLOBAL::ProblemType::gas_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::PartFS3I1Wc(comm));
    }
    break;
    case GLOBAL::ProblemType::ac_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::ACFSI(comm));
    }
    break;
    case GLOBAL::ProblemType::thermo_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::PartFS3I2Wc(comm));
    }
    break;
    case GLOBAL::ProblemType::biofilm_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::BiofilmFSI(comm));
    }
    break;
    case GLOBAL::ProblemType::fps3i:
    {
      fs3i = Teuchos::rcp(new FS3I::PartFpS3I1Wc(comm));
    }
    break;
    default:
      FOUR_C_THROW("solution of unknown problemtyp %d requested", probtype);
      break;
  }

  fs3i->Init();
  fs3i->Setup();

  // read the restart information, set vectors and variables ---
  // be careful, dofmaps might be changed here in a Redistribute call
  fs3i->ReadRestart();

  // if running FPS3I in parallel one needs to redistribute the interface after restarting
  fs3i->redistribute_interface();

  // now do the coupling and create combined dofmaps
  fs3i->SetupSystem();

  fs3i->Timeloop();

  fs3i->TestResults(comm);

  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = CORE::COMM::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
}

FOUR_C_NAMESPACE_CLOSE
