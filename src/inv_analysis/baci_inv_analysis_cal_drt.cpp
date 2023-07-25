/*----------------------------------------------------------------------*/
/*! \file
\brief Start routines for the inverse analysis

\level 3

*/
/*----------------------------------------------------------------------*/
#include "baci_lib_globalproblem.H"
#include "baci_inv_analysis_cal_drt.H"
#include "baci_inpar_statinvanalysis.H"
#include "baci_inpar_invanalysis.H"

// the new inverse analysis
#include "baci_inv_analysis_control.H"
#include "baci_inv_analysis_base.H"

/*----------------------------------------------------------------------*/
void invana_cal()
{
  // get input lists
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  if (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvAnalysisType>(invp, "STAT_INV_ANALYSIS") !=
      INPAR::INVANA::stat_inv_none)
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
