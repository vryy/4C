/*---------------------------------------------------------------------*/
/*! \file

\brief A map extractor for Dirichlet conditions

\level 0


*/
/*---------------------------------------------------------------------*/

#include "lib_dirichletextractor.H"

#include "lib_discret.H"
#include "lib_condition_selector.H"
#include "utils_exceptions.H"


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
