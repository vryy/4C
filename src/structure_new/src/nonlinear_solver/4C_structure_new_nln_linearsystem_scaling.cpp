// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_nln_linearsystem_scaling.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_structure_new_timint_basedatasdyn.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Teuchos_ParameterList.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Solid::Nln::LinSystem::StcScaling::StcScaling(
    const Solid::TimeInt::BaseDataSDyn& DataSDyn, Solid::TimeInt::BaseDataGlobalState& GState)
    : stcscale_(DataSDyn.get_stc_algo_type()), stclayer_(DataSDyn.get_stc_layer()), stcmat_(nullptr)
{
  // prepare matrix for scaled thickness business of thin shell structures
  stcmat_ =
      std::make_shared<Core::LinAlg::SparseMatrix>(*GState.dof_row_map_view(), 81, true, true);
  stcmat_->zero();

  // create the parameters for the discretization
  Teuchos::ParameterList p;

  // get discretization
  std::shared_ptr<Core::FE::Discretization> discret = GState.get_discret();

  // action for elements
  discret->set_state("displacement", GState.get_dis_np());

  const std::string action = "calc_stc_matrix";
  p.set("action", action);
  p.set<Inpar::Solid::StcScale>("stc_scaling", stcscale_);
  p.set("stc_layer", 1);

  discret->evaluate(p, stcmat_, nullptr, nullptr, nullptr, nullptr);

  stcmat_->complete();

  for (int lay = 2; lay <= stclayer_; ++lay)
  {
    Teuchos::ParameterList pe;

    pe.set("action", action);
    pe.set<Inpar::Solid::StcScale>("stc_scaling", stcscale_);
    pe.set("stc_layer", lay);

    std::shared_ptr<Core::LinAlg::SparseMatrix> tmpstcmat =
        std::make_shared<Core::LinAlg::SparseMatrix>(*GState.dof_row_map_view(), 81, true, true);
    tmpstcmat->zero();

    discret->evaluate(pe, tmpstcmat, nullptr, nullptr, nullptr, nullptr);
    tmpstcmat->complete();

    stcmat_ = Core::LinAlg::matrix_multiply(*tmpstcmat, false, *stcmat_, false, true, false, true);
  }

  discret->clear_state();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::Nln::LinSystem::StcScaling::scaleLinearSystem(Epetra_LinearProblem& problem)
{
  // get stiffness matrix
  Epetra_CrsMatrix* stiffmat = dynamic_cast<Epetra_CrsMatrix*>(problem.GetMatrix());
  std::shared_ptr<Epetra_CrsMatrix> stiff_epetra = Core::Utils::shared_ptr_from_ref(*stiffmat);
  std::shared_ptr<Core::LinAlg::SparseMatrix> stiff_linalg =
      std::make_shared<Core::LinAlg::SparseMatrix>(stiff_epetra, Core::LinAlg::View);

  // get rhs
  Core::LinAlg::VectorView rhs_view(*dynamic_cast<Epetra_Vector*>(problem.GetRHS()));
  Core::LinAlg::Vector<double>& rhs(rhs_view);

  // right multiplication of stiffness matrix
  stiff_scaled_ =
      Core::LinAlg::matrix_multiply(*stiff_linalg, false, *stcmat_, false, true, false, true);

  // left multiplication of stiffness matrix and rhs
  if (stcscale_ == Inpar::Solid::stc_currsym)
  {
    stiff_scaled_ =
        Core::LinAlg::matrix_multiply(*stcmat_, true, *stiff_scaled_, false, true, false, true);

    std::shared_ptr<Core::LinAlg::Vector<double>> rhs_scaled =
        Core::LinAlg::create_vector(problem.GetRHS()->Map(), true);
    stcmat_->multiply(true, rhs, *rhs_scaled);
    rhs.Update(1.0, *rhs_scaled, 0.0);
  }

  // set new stiffness matrix
  problem.SetOperator(stiff_scaled_->epetra_matrix().get());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::Nln::LinSystem::StcScaling::unscaleLinearSystem(Epetra_LinearProblem& problem)
{
  std::shared_ptr<Core::LinAlg::Vector<double>> disisdc =
      Core::LinAlg::create_vector(problem.GetLHS()->Map(), true);
  Epetra_MultiVector* disi = problem.GetLHS();

  Core::LinAlg::VectorView disi_view(*disi);
  stcmat_->multiply(false, disi_view, *disisdc);
  disi->Update(1.0, *disisdc, 0.0);
}

FOUR_C_NAMESPACE_CLOSE
