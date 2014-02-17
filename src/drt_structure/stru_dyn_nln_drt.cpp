/*!----------------------------------------------------------------------

\file stru_dyn_nln_drt.cpp

\brief Control routine for structural dynamics (outsourced to adapter layer)

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/ad_str_structure.H"
#include "stru_dyn_nln_drt.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "stru_resulttest.H"
#include "str_invanalysis.H"
#include "str_statinvanalysis.H"
#include "../drt_inv_analysis/inv_analysis.H"
#include "../drt_inv_analysis/stat_inv_analysis.H"

#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

//periodic boundary conditions
#include "../drt_fluid/drt_periodicbc.H"

#ifdef HAVE_FFTW
#include "str_mlmc.H"
#endif

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void caldyn_drt()
{
  // get input lists
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  //get list for multi level monte carlo
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // do we want to do inverse analysis?
  if (DRT::INPUT::IntegralValue<INPAR::STR::InvAnalysisType>(iap,"INV_ANALYSIS")
      != INPAR::STR::inv_none)
  {
    STR::invanalysis();
  }
  // do we want to do statistical inverse analysis?
  else if (DRT::INPUT::IntegralValue<INPAR::STR::StatInvAnalysisType>(statinvp,"STAT_INV_ANALYSIS")
     != INPAR::STR::stat_inv_none)
  {
    STR::statinvanalysis();
  }
  //do we want multi level monte carlo
  else if (Teuchos::getIntegralValue<int>(mlmcp,"MLMC")!= false)
  {
#ifdef HAVE_FFTW
    STR::mlmc();
#else
  std::cout<< RED_LIGHT << "CANNOT PERFORM MLMC WITHOUT FFTW  "<< END_COLOR << std::endl;
#endif
  }
  else
  {
    // get input lists
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    // major switch to different time integrators
    switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
    {
    case INPAR::STR::dyna_statics:
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_onesteptheta:
    case INPAR::STR::dyna_gemm:
    case INPAR::STR::dyna_expleuler:
    case INPAR::STR::dyna_centrdiff:
    case INPAR::STR::dyna_ab2:
    case INPAR::STR::dyna_euma:
    case INPAR::STR::dyna_euimsto:
    case INPAR::STR::dyna_statmech:
      dyn_nlnstructural_drt();
      break;
    default:
      dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
      break;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | structural nonlinear dynamics                                        |
 *----------------------------------------------------------------------*/
void dyn_nlnstructural_drt()
{
  // get input lists
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

  // connect degrees of freedom for periodic boundary conditions
  {
    PeriodicBoundaryConditions pbc_struct(structdis);

    if (pbc_struct.HasPBC())
    {
      pbc_struct.UpdateDofsForPeriodicBoundaryConditions();
    }
  }

  // create an adapterbase and adapter
  ADAPTER::StructureBaseAlgorithm adapterbase(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
  ADAPTER::Structure& structadaptor = const_cast<ADAPTER::Structure&>(adapterbase.StructureField());

  // do restart
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    structadaptor.ReadRestart(restart);
  }

  // write output at beginnning of calc
  else
  {
    //Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->GetDis("structure");
    //Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
    //output->NewStep(0, 0.0);
    //Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(*(actdis->DofRowMap())));
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
  Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(structdis->Comm());
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, true);

  // time to go home...
  return;

} // end of dyn_nlnstructural_drt()

