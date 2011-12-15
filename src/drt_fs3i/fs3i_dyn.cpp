
#ifdef CCADISCRET

#include "fs3i.H"
#include "gas_fsi.H"
#include "biofilm_fsi.H"
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
  const Epetra_MpiComm& epetrampicomm = dynamic_cast<const Epetra_MpiComm&>(DRT::Problem::Instance()->Dis(genprob.numff,0)->Comm());
  if (!(&epetrampicomm))
    dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  Epetra_MpiComm& comm = const_cast<Epetra_MpiComm&>(epetrampicomm);
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
