
#ifdef CCADISCRET

#include "fs3i.H"
#include "gas_fsi.H"
#include "biofilm_fsi.H"
#include "thermo_fsi.H"
#include "aero_tfsi.H"
#include "fs3i_dyn.H"
#include "../drt_lib/drt_globalproblem.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
// entry point for all kinds if FS3I
/*----------------------------------------------------------------------*/
void fs3i_dyn()
{
#ifdef PARALLEL
  const Epetra_Comm& comm = DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm();
#else
  Epetra_SerialComm comm;
#endif

  Teuchos::RCP<FS3I::FS3I_Base> fs3i;

  switch (genprob.probtyp)
  {
    case prb_gas_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::GasFSI(comm));
    }
    break;
    case prb_biofilm_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::BiofilmFSI(comm));
    }
    break;
    case prb_thermo_fsi:
    {
      fs3i = Teuchos::rcp(new FS3I::ThermoFSI(comm));
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

  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
}

#endif
