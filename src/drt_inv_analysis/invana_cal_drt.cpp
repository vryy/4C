/*----------------------------------------------------------------------*/
/*! \file
\brief Start routines for the inverse analysis

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "invana_cal_drt.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_inpar/inpar_invanalysis.H"

// the hook to the old inverse analyses
#include "../drt_structure/str_invanalysis.H"

// the new inverse analysis
#include "invana_control.H"
#include "invana_base.H"

/*----------------------------------------------------------------------*/
void invana_cal()
{
  // get input lists
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();

  if ((DRT::INPUT::IntegralValue<INPAR::STR::InvAnalysisType>(iap, "INV_ANALYSIS") !=
          INPAR::STR::inv_none) and
      (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvAnalysisType>(invp, "STAT_INV_ANALYSIS") !=
          INPAR::INVANA::stat_inv_none))
    dserror("Use either INV_ANALYSIS or STAT_INV_ANALYSIS. They cannot be used simultaneously.");
  // do we want to do inverse analysis?
  else if (DRT::INPUT::IntegralValue<INPAR::STR::InvAnalysisType>(iap, "INV_ANALYSIS") !=
           INPAR::STR::inv_none)
  {
    STR::invanalysis();
  }
  else if (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvAnalysisType>(
               invp, "STAT_INV_ANALYSIS") != INPAR::INVANA::stat_inv_none)
  {
    // initialize inverse solution process
    INVANA::InvanaControl inversesolution;
    inversesolution.Init(invp);

    // solve
    int restart = DRT::Problem::Instance()->Restart();
    inversesolution.Solve(restart);

    // test
    DRT::Problem::Instance()->AddFieldTest(inversesolution.CreateFieldTest());
    DRT::Problem::Instance()->TestAll(inversesolution.InverseProblem()->Comm());
  }

  return;
}
