/*!----------------------------------------------------------------------

\file invana_cal_drt.H

\brief Start routines for the inverse analysis

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10361
</pre>

*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "invana_cal_drt.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inv_analysis/invana_factory.H"
#include "../drt_inv_analysis/optimizer_factory.H"
#include "../drt_inv_analysis/invana_base.H"
#include "../drt_structure/str_invanalysis.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void invana_cal()
{
  // get input lists
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();

  // do we want to do inverse analysis?
  if (DRT::INPUT::IntegralValue<INPAR::STR::InvAnalysisType>(iap,"INV_ANALYSIS")
      != INPAR::STR::inv_none)
  {
    STR::invanalysis();
  }
  else if (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvAnalysisType>(invp,"STAT_INV_ANALYSIS")
   != INPAR::INVANA::stat_inv_none)
  {

    // access the discretization
    Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
    actdis = DRT::Problem::Instance()->GetDis("structure");

    // set degrees of freedom in the discretization
    if (!actdis->Filled()) actdis->FillComplete();
    if (!actdis->HaveDofs()) actdis->FillComplete();

    // create an instance of an optimization problem
    INVANA::InvanaFactory invfac;
    Teuchos::RCP<INVANA::InvanaBase> optprob = invfac.Create(actdis,invp);

    // solve
    int restart= DRT::Problem::Instance()->Restart();
    optprob->Solve(restart);

    // test
    DRT::Problem::Instance()->AddFieldTest(optprob->CreateFieldTest());
    DRT::Problem::Instance()->TestAll(actdis->Comm());
  }

  return;
}
