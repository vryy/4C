/*---------------------------------------------------------------------*/
/*!

\brief Utils methods concerning the discretization evaluation

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "drt_utils_discret.H"
#include "drt_assemblestrategy.H"
#include "drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Evaluate(DRT::Discretization& discret, Teuchos::ParameterList& eparams,
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    const Teuchos::RCP<Epetra_Vector>& systemvector, const Epetra_Map* col_ele_map)
{
  std::vector<Teuchos::RCP<LINALG::SparseOperator>> systemmatrices(2, Teuchos::null);
  std::vector<Teuchos::RCP<Epetra_Vector>> systemvectors(3, Teuchos::null);

  systemmatrices[0] = systemmatrix;
  systemvectors[0] = systemvector;

  Evaluate(discret, eparams, systemmatrices, systemvectors, col_ele_map);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Evaluate(DRT::Discretization& discret, Teuchos::ParameterList& eparams,
    std::vector<Teuchos::RCP<LINALG::SparseOperator>>& systemmatrices,
    std::vector<Teuchos::RCP<Epetra_Vector>>& systemvectors, const Epetra_Map* col_ele_map)
{
  dsassert(systemmatrices.size() <= 2,
      "Currently a maximum number of two "
      "system-matrices is supported!");
  dsassert(systemvectors.size() <= 3,
      "Currently a maximum number of three "
      "system-vectors is supported!");

  if (systemmatrices.size() < 2) systemmatrices.resize(2, Teuchos::null);
  if (systemvectors.size() < 3) systemvectors.resize(3, Teuchos::null);

  DRT::AssembleStrategy strategy(0, 0, systemmatrices[0], systemmatrices[1], systemvectors[0],
      systemvectors[1], systemvectors[2]);
  Evaluate(discret, eparams, strategy, col_ele_map);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Evaluate(DRT::Discretization& discret, Teuchos::ParameterList& eparams,
    DRT::AssembleStrategy& strategy, const Epetra_Map* col_ele_map)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::UTILS::Evaluate");

  if (!discret.Filled()) dserror("FillComplete() was not called");
  if (!discret.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  int row = strategy.FirstDofSet();
  int col = strategy.SecondDofSet();

  // call the element's register class pre-evaluation method
  // for each type of element
  // for most element types, just the base class dummy is called
  // that does nothing
  {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::UTILS::Evaluate PreEvaluate");
    ParObjectFactory::Instance().PreEvaluate(discret, eparams, strategy.Systemmatrix1(),
        strategy.Systemmatrix2(), strategy.Systemvector1(), strategy.Systemvector2(),
        strategy.Systemvector3());
  }

  Element::LocationArray la(discret.NumDofSets());

  bool is_subset = false;
  if (not col_ele_map)
    col_ele_map = discret.ElementColMap();
  else
    is_subset = true;

  // loop over column elements
  const int numcolele = col_ele_map->NumMyElements();
  const int* ele_gids = col_ele_map->MyGlobalElements();

  for (int i = 0; i < numcolele; ++i)
  {
    DRT::Element* actele = NULL;
    if (is_subset)
    {
      const int egid = ele_gids[i];
      actele = discret.gElement(egid);
    }
    else
      actele = discret.lColElement(i);

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::UTILS::Evaluate LocationVector");
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(discret, la, false);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::UTILS::Evaluate Resize");

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage(la[row].Size(), la[col].Size());
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::UTILS::Evaluate elements");
      // call the element evaluate method
      int err = actele->Evaluate(eparams, discret, la, strategy.Elematrix1(), strategy.Elematrix2(),
          strategy.Elevector1(), strategy.Elevector2(), strategy.Elevector3());
      if (err)
        dserror("Proc %d: Element %d returned err=%d", discret.Comm().MyPID(), actele->Id(), err);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::UTILS::Evaluate assemble");
      int eid = actele->Id();
      strategy.AssembleMatrix1(eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_);
      strategy.AssembleMatrix2(eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_);
      strategy.AssembleVector1(la[row].lm_, la[row].lmowner_);
      strategy.AssembleVector2(la[row].lm_, la[row].lmowner_);
      strategy.AssembleVector3(la[row].lm_, la[row].lmowner_);
    }

  }  // loop over all considered elements

  return;
}
