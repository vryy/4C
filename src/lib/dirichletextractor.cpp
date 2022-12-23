/*---------------------------------------------------------------------*/
/*! \file

\brief A map extractor for Dirichlet conditions

\level 0


*/
/*---------------------------------------------------------------------*/

#include "dirichletextractor.H"

#include "discret.H"
#include "condition_selector.H"
#include "condition_utils.H"
#include "dserror.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::DirichletExtractor::Setup(const DRT::Discretization& dis)
{
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::DirichletSelector(dis)));
  mcs.SetupExtractor(dis, *dis.DofRowMap(), *this);
}

void DRT::DirichletExtractor::ZeroDirichlets(Teuchos::RCP<Epetra_Vector> residual) const
{
  DirichletPutScalar(*residual, 0.0);
}
