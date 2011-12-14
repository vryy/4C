/*!----------------------------------------------------------------------

\file stru_dyn_nln_drt.cpp

\brief Control routine for structural dynamics (outsourced to adapter layer)

<pre>
Maintainer: Thomas Klöppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_dyn_nln_drt.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "stru_resulttest.H"
#include "str_invanalysis.H"
#include "../drt_inv_analysis/inv_analysis.H"

#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

#ifdef HAVE_FFTW
#include "str_mlmc.H"
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
extern "C"
void caldyn_drt()
{
  // get input lists
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  //get list for multi level monte carlo
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // do we want to do inverse analysis?
  if (DRT::INPUT::IntegralValue<INPAR::STR::InvAnalysisType>(iap,"INV_ANALYSIS")
      != INPAR::STR::inv_none)
  {
    STR::invanalysis();
  }
  //
  //do we want multi level monte carlo
  else if (Teuchos::getIntegralValue<int>(mlmcp,"MLMC")!= false)
  {
#ifdef HAVE_FFTW
    STR::mlmc();
#else
  cout<< RED_LIGHT << "CANNOT PERFORM MLMC WITHOUT FFTW  "<< END_COLOR << endl;
#endif
  }
  else
  {
    // get input lists
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

    // major switch to different time integrators
    switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
    {
    case INPAR::STR::dyna_centr_diff:
      dserror("no central differences in DRT");
      break;
    case INPAR::STR::dyna_gen_alfa:
    case INPAR::STR::dyna_gen_alfa_statics:
    case INPAR::STR::dyna_statics:
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_onesteptheta:
    case INPAR::STR::dyna_gemm:
    case INPAR::STR::dyna_ab2:
    case INPAR::STR::dyna_euma :
    case INPAR::STR::dyna_euimsto :
      dyn_nlnstructural_drt();
      break;
    case INPAR::STR::dyna_Gen_EMM:
      dserror("GEMM not supported");
      break;
    default:
      dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
    }
  }
}


/*----------------------------------------------------------------------*
 | structural nonlinear dynamics                                        |
 *----------------------------------------------------------------------*/
void dyn_nlnstructural_drt()
{
  // create an adapterbase and adapter
  ADAPTER::StructureBaseAlgorithm adapterbase(DRT::Problem::Instance()->StructuralDynamicParams());
  ADAPTER::Structure& structadaptor = const_cast<ADAPTER::Structure&>(adapterbase.StructureField());

  // do restart
  if (genprob.restart)
  {
    structadaptor.ReadRestart(genprob.restart);
  }

  // write output at beginnning of calc
  else
  {
    //RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);
    //RCP<IO::DiscretizationWriter> output = rcp(new IO::DiscretizationWriter(actdis));
    //output->NewStep(0, 0.0);
    //RCP<Epetra_Vector> zeros = rcp (new Epetra_Vector(*(actdis->DofRowMap())));
    //output->WriteVector("displacement",zeros);
    //output->WriteElementData();
  }

#if 1
  // run time integration
  structadaptor.Integrate();
#else // this is a bone optimization hack - do not touch
  structadaptor.PrepareTimeStep();
  structadaptor.Solve();
  structadaptor.Update();
  structadaptor.Output();

  const Epetra_Vector* disp = structadaptor.Dispn().get();
  const_cast<Epetra_Vector*>(disp)->PutScalar(0.0);

  structadaptor.PrepareTimeStep();
  structadaptor.Solve();
  structadaptor.Update();
  structadaptor.Output();
#endif

  // test results
  DRT::Problem::Instance()->AddFieldTest(structadaptor.CreateFieldTest());
  DRT::Problem::Instance()->TestAll(structadaptor.DofRowMap()->Comm());

  // print monitoring of time consumption
  Teuchos::TimeMonitor::summarize();

  // time to go home...
  return;

} // end of dyn_nlnstructural_drt()

#endif  // #ifdef CCADISCRET
