
#ifdef CCADISCRET

#include "fs3i_dyn.H"

#include "fs3i.H"
#include "fs3i_partitioned_1wc.H"
#include "fs3i_partitioned_2wc.H"
#include "biofilm_fsi.H"
#include "aero_tfsi.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_comm/comm_utils.H"

#include <Teuchos_TimeMonitor.hpp>

extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
// entry point for all kinds of FS3I
/*----------------------------------------------------------------------*/
void fs3i_dyn()
{
  const Epetra_Comm& comm = DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm();

  Teuchos::RCP<FS3I::FS3I_Base> fs3i;

  switch (genprob.probtyp)
  {
    case prb_gas_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::PartFS3I_1WC(comm));
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
    case prb_tfsi_aero:
    {
      fs3i = Teuchos::rcp(new FS3I::AeroTFSI(comm));
    }
    break;
    default:
      dserror("solution of unknown problemtyp %d requested", genprob.probtyp);
    break;
  }

  // read the restart information, set vectors and variables ---
  // be careful, dofmaps might be changed here in a Redistribute call
  fs3i->ReadRestart();

  // now do the coupling and create combined dofmaps
  fs3i->SetupSystem();

  fs3i->Timeloop();

  fs3i->TestResults(comm);

  Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
}

#endif
