/*----------------------------------------------------------------------*/
/*!
\file str_statinvanalysis.cpp
\brief Statistical inverse analysis for structures

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de

</pre>
*/

#include "../drt_lib/drt_globalproblem.H"
#include "str_statinvanalysis.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inv_analysis/invana_factory.H"
#include "../drt_inv_analysis/optimizer_factory.H"
#include "../drt_inv_analysis/invana_base.H"

/*======================================================================*/
/* inverse analysis of structures */
void STR::statinvanalysis()
{
  // get input lists
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("structure");

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();
  if (!actdis->HaveDofs()) actdis->FillComplete();

  // create an instance of an optimization problem
  STR::INVANA::InvanaFactory invfac;
  Teuchos::RCP<STR::INVANA::InvanaBase> optprob = invfac.Create(actdis,invp);

  // solve
  int restart= DRT::Problem::Instance()->Restart();
  optprob->Solve(restart);

  // test
  DRT::Problem::Instance()->AddFieldTest(optprob->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(actdis->Comm());

  // done
  return;
} // end str_statinvanalysis()


/*----------------------------------------------------------------------*/
