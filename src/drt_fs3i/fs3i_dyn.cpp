
#ifdef CCADISCRET

#include "fs3i.H"
#include "fs3i_1wc.H"
#include "fs3i_biofilm_growth.H"
#include "fs3i_dyn.H"
#include "../drt_lib/drt_globalproblem.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
// entry point for all kinds if FS3I
/*----------------------------------------------------------------------*/
void fs3i_dyn()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Teuchos::RCP<FS3I::FS3I_Base> fs3i;

  switch (genprob.probtyp)
  {
    case prb_fsi_lung_gas:
    {
      // this has to be changed! -> introduce new problem type for biofilm!
      const Teuchos::ParameterList& biofilmcontrol = DRT::Problem::Instance()->BIOFILMControlParams();
      const int surfgrowth = DRT::INPUT::IntegralValue<int>(biofilmcontrol,"SURFACEGROWTH");

      if (!surfgrowth)
        fs3i = Teuchos::rcp(new FS3I::FS3I_1WC(comm));
      else
      {
        fs3i = Teuchos::rcp(new FS3I::BiofilmGrowth(comm));
      }
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
