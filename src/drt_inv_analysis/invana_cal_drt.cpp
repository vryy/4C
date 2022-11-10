/*----------------------------------------------------------------------*/
/*! \file
\brief Start routines for the inverse analysis

\level 3

*/
/*----------------------------------------------------------------------*/
#include "drt_globalproblem.H"
#include "invana_cal_drt.H"
#include "inpar_statinvanalysis.H"
#include "inpar_invanalysis.H"

// the new inverse analysis
#include "invana_control.H"
#include "invana_base.H"

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
