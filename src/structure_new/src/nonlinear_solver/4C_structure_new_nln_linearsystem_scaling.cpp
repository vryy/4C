/*-----------------------------------------------------------*/
/*! \file

\brief


\level 3

 */
/*-----------------------------------------------------------*/

#include "4C_structure_new_nln_linearsystem_scaling.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_structure_new_timint_basedatasdyn.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Solid::Nln::LinSystem::StcScaling::StcScaling(
    const Solid::TimeInt::BaseDataSDyn& DataSDyn, Solid::TimeInt::BaseDataGlobalState& GState)
    : stcscale_(DataSDyn.get_stc_algo_type()),
      stclayer_(DataSDyn.get_stc_layer()),
      stcmat_(Teuchos::null)
{
  // prepare matrix for scaled thickness business of thin shell structures
  stcmat_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*GState.dof_row_map_view(), 81, true, true));
  stcmat_->zero();

  // create the parameters for the discretization
  Teuchos::ParameterList p;

  // get discretization
  Teuchos::RCP<Core::FE::Discretization> discret = GState.get_discret();

  // action for elements
  discret->set_state("displacement", GState.get_dis_np());

  const std::string action = "calc_stc_matrix";
  p.set("action", action);
  p.set<int>("stc_scaling", stcscale_);
  p.set("stc_layer", 1);

  discret->evaluate(p, stcmat_, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  stcmat_->complete();

  for (int lay = 2; lay <= stclayer_; ++lay)
  {
    Teuchos::ParameterList pe;

    pe.set("action", action);
    pe.set<int>("stc_scaling", stcscale_);
    pe.set("stc_layer", lay);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> tmpstcmat =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*GState.dof_row_map_view(), 81, true, true));
    tmpstcmat->zero();

    discret->evaluate(pe, tmpstcmat, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    tmpstcmat->complete();

    stcmat_ = MLMultiply(*tmpstcmat, *stcmat_, true, false, true);
  }

  discret->clear_state();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::Nln::LinSystem::StcScaling::scaleLinearSystem(Epetra_LinearProblem& problem)
{
  // get stiffness matrix
  Epetra_CrsMatrix* stiffmat = dynamic_cast<Epetra_CrsMatrix*>(problem.GetMatrix());
  Teuchos::RCP<Epetra_CrsMatrix> stiff_epetra = Teuchos::rcp(stiffmat, false);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_linalg =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(stiff_epetra, Core::LinAlg::View));

  // get rhs
  Epetra_Vector* rhs = dynamic_cast<Epetra_Vector*>(problem.GetRHS());

  // right multiplication of stiffness matrix
  stiff_scaled_ = MLMultiply(*stiff_linalg, *stcmat_, true, false, true);

  // left multiplication of stiffness matrix and rhs
  if (stcscale_ == Inpar::Solid::stc_currsym)
  {
    stiff_scaled_ = MLMultiply(*stcmat_, true, *stiff_scaled_, false, true, false, true);

    Teuchos::RCP<Epetra_Vector> rhs_scaled =
        Core::LinAlg::CreateVector(problem.GetRHS()->Map(), true);
    stcmat_->multiply(true, *rhs, *rhs_scaled);
    rhs->Update(1.0, *rhs_scaled, 0.0);
  }

  // set new stiffness matrix
  problem.SetOperator(stiff_scaled_->epetra_matrix().get());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::Nln::LinSystem::StcScaling::unscaleLinearSystem(Epetra_LinearProblem& problem)
{
  Teuchos::RCP<Epetra_MultiVector> disisdc =
      Core::LinAlg::CreateVector(problem.GetLHS()->Map(), true);
  Epetra_MultiVector* disi = dynamic_cast<Epetra_Vector*>(problem.GetLHS());

  stcmat_->multiply(false, *disi, *disisdc);
  disi->Update(1.0, *disisdc, 0.0);
}

FOUR_C_NAMESPACE_CLOSE
