/*-----------------------------------------------------------*/
/*!
\file str_nln_linearsystem_scaling.cpp

\brief

\maintainer Amadeus Gebauer

\level 3

 */
/*-----------------------------------------------------------*/

#include "str_nln_linearsystem_scaling.H"

#include <iostream>
#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>

#include "../drt_lib/drt_discret_interface.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_structure_new/str_timint_basedatasdyn.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::NLN::LinSystem::StcScaling::StcScaling(
    STR::TIMINT::BaseDataSDyn DataSDyn, STR::TIMINT::BaseDataGlobalState GState)
    : stcscale_(DataSDyn.GetSTCAlgoType()),
      stclayer_(DataSDyn.GetSTCLayer()),
      stcmat_(Teuchos::null)
{
  // prepare matrix for scaled thickness business of thin shell structures
  stcmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*GState.DofRowMapView(), 81, true, true));
  stcmat_->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList p;

  // get discretization
  Teuchos::RCP<DRT::DiscretizationInterface> discret = GState.GetMutableDiscret();

  // action for elements
  discret->SetState("displacement", GState.GetDisNp());

  const std::string action = "calc_stc_matrix";
  p.set("action", action);
  p.set<int>("stc_scaling", stcscale_);
  p.set("stc_layer", 1);

  discret->Evaluate(p, stcmat_, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  stcmat_->Complete();

  for (int lay = 2; lay <= stclayer_; ++lay)
  {
    Teuchos::ParameterList pe;

    pe.set("action", action);
    pe.set<int>("stc_scaling", stcscale_);
    pe.set("stc_layer", lay);

    Teuchos::RCP<LINALG::SparseMatrix> tmpstcmat =
        Teuchos::rcp(new LINALG::SparseMatrix(*GState.DofRowMapView(), 81, true, true));
    tmpstcmat->Zero();

    discret->Evaluate(pe, tmpstcmat, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    tmpstcmat->Complete();

    stcmat_ = MLMultiply(*tmpstcmat, *stcmat_, true, false, true);
  }

  discret->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::NLN::LinSystem::StcScaling::scaleLinearSystem(Epetra_LinearProblem& problem)
{
  // get stiffness matrix
  Epetra_CrsMatrix* stiffmat = dynamic_cast<Epetra_CrsMatrix*>(problem.GetMatrix());
  Teuchos::RCP<Epetra_CrsMatrix> stiff_epetra = Teuchos::rcp(stiffmat, false);
  Teuchos::RCP<LINALG::SparseMatrix> stiff_linalg =
      Teuchos::rcp(new LINALG::SparseMatrix(stiff_epetra, LINALG::View));

  // get rhs
  Epetra_Vector* rhs = dynamic_cast<Epetra_Vector*>(problem.GetRHS());

  // right multiplication of stiffness matrix
  stiff_scaled_ = MLMultiply(*stiff_linalg, *stcmat_, true, false, true);

  // left multiplication of stiffness matrix and rhs
  if (stcscale_ == INPAR::STR::stc_currsym)
  {
    stiff_scaled_ = MLMultiply(*stcmat_, true, *stiff_scaled_, false, true, false, true);

    Teuchos::RCP<Epetra_Vector> rhs_scaled = LINALG::CreateVector(problem.GetRHS()->Map(), true);
    stcmat_->Multiply(true, *rhs, *rhs_scaled);
    rhs->Update(1.0, *rhs_scaled, 0.0);
  }

  // set new stiffness matrix
  problem.SetOperator(stiff_scaled_->EpetraMatrix().get());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::NLN::LinSystem::StcScaling::unscaleLinearSystem(Epetra_LinearProblem& problem)
{
  Teuchos::RCP<Epetra_MultiVector> disisdc = LINALG::CreateVector(problem.GetLHS()->Map(), true);
  Epetra_MultiVector* disi = dynamic_cast<Epetra_Vector*>(problem.GetLHS());

  stcmat_->Multiply(false, *disi, *disisdc);
  disi->Update(1.0, *disisdc, 0.0);
}
