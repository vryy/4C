/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche poro contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_nitsche_strategy_poro.hpp"

#include "4C_contact_interface.hpp"
#include "4C_contact_nitsche_utils.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_so3_plast_ssn.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_Operator.h>

FOUR_C_NAMESPACE_OPEN

void CONTACT::NitscheStrategyPoro::apply_force_stiff_cmt(
    Teuchos::RCP<Core::LinAlg::Vector<double>> dis, Teuchos::RCP<Core::LinAlg::SparseOperator>& kt,
    Teuchos::RCP<Core::LinAlg::Vector<double>>& f, const int step, const int iter, bool predictor)
{
  if (predictor) return;

  CONTACT::NitscheStrategy::apply_force_stiff_cmt(dis, kt, f, step, iter, predictor);

  // Evaluation for all interfaces
  fp_ = create_rhs_block_ptr(CONTACT::VecBlockType::porofluid);
  kpp_ = create_matrix_block_ptr(CONTACT::MatBlockType::porofluid_porofluid);
  kpd_ = create_matrix_block_ptr(CONTACT::MatBlockType::porofluid_displ);
  kdp_ = create_matrix_block_ptr(CONTACT::MatBlockType::displ_porofluid);
  //    for (int i = 0; i < (int) interface_.size(); ++i)
  //    {
  //      for (int e=0;e<interface_[i]->discret().ElementColMap()->NumMyElements();++e)
  //      {
  //        Mortar::Element* mele
  //        =dynamic_cast<Mortar::Element*>(interface_[i]->discret().gElement(
  //            interface_[i]->discret().ElementColMap()->GID(e)));
  //        mele->get_nitsche_container().ClearAll();
  //      }
  //    }
}

void CONTACT::NitscheStrategyPoro::set_state(
    const enum Mortar::StateType& statename, const Core::LinAlg::Vector<double>& vec)
{
  if (statename == Mortar::state_svelocity)
  {
    Teuchos::RCP<Core::FE::Discretization> dis = Global::Problem::instance()->get_dis("structure");
    if (dis == Teuchos::null) FOUR_C_THROW("didn't get my discretization");
    set_parent_state(statename, vec, *dis);
  }
  else
    CONTACT::NitscheStrategy::set_state(statename, vec);
}

void CONTACT::NitscheStrategyPoro::set_parent_state(const enum Mortar::StateType& statename,
    const Core::LinAlg::Vector<double>& vec, const Core::FE::Discretization& dis)
{
  //
  if (statename == Mortar::state_fvelocity || statename == Mortar::state_fpressure)
  {
    Core::LinAlg::Vector<double> global(*dis.dof_col_map(), true);
    Core::LinAlg::export_to(vec, global);

    // set state on interfaces
    for (const auto& interface : interface_)
    {
      Core::FE::Discretization& idiscret = interface->discret();

      for (int j = 0; j < interface->discret().element_col_map()->NumMyElements(); ++j)
      {
        const int gid = interface->discret().element_col_map()->GID(j);

        auto* ele = dynamic_cast<Mortar::Element*>(idiscret.g_element(gid));

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;

        if (ele->parent_slave_element())  // if this pointer is nullptr, this parent is impermeable
        {
          // this gets values in local order
          ele->parent_slave_element()->location_vector(dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          Core::FE::extract_my_values(global, myval, lm);

          std::vector<double> vel;
          std::vector<double> pres;

          for (int n = 0; n < ele->parent_slave_element()->num_node(); ++n)
          {
            for (unsigned dim = 0; dim < 3; ++dim)
            {
              vel.push_back(myval[n * 4 + dim]);
            }
            pres.push_back(myval[n * 4 + 3]);
          }

          ele->mo_data().parent_pf_pres() = pres;
          ele->mo_data().parent_pf_vel() = vel;
          ele->mo_data().parent_pf_dof() = lm;
        }
      }
    }
  }
  else
    CONTACT::NitscheStrategy::set_parent_state(statename, vec, dis);
}

Teuchos::RCP<Epetra_FEVector> CONTACT::NitscheStrategyPoro::setup_rhs_block_vec(
    const enum CONTACT::VecBlockType& bt) const
{
  switch (bt)
  {
    case CONTACT::VecBlockType::porofluid:
      return Teuchos::make_rcp<Epetra_FEVector>(
          *Global::Problem::instance()->get_dis("porofluid")->dof_row_map());
    default:
      return CONTACT::NitscheStrategy::setup_rhs_block_vec(bt);
  }
}

Teuchos::RCP<const Core::LinAlg::Vector<double>> CONTACT::NitscheStrategyPoro::get_rhs_block_ptr(
    const enum CONTACT::VecBlockType& bp) const
{
  if (!curr_state_eval_) FOUR_C_THROW("you didn't evaluate this contact state first");

  switch (bp)
  {
    case CONTACT::VecBlockType::porofluid:

      return Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*fp_);
    default:
      return CONTACT::NitscheStrategy::get_rhs_block_ptr(bp);
  }
}

Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategyPoro::setup_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_porofluid:
      return Teuchos::make_rcp<Core::LinAlg::SparseMatrix>(
          *Global::Problem::instance()->get_dis("structure")->dof_row_map(), 100, true, false,
          Core::LinAlg::SparseMatrix::FE_MATRIX);
    case CONTACT::MatBlockType::porofluid_displ:
    case CONTACT::MatBlockType::porofluid_porofluid:
      return Teuchos::make_rcp<Core::LinAlg::SparseMatrix>(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *Global::Problem::instance()->get_dis("porofluid")->dof_row_map()),
          100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    default:
      return CONTACT::NitscheStrategy::setup_matrix_block_ptr(bt);
  }
}

void CONTACT::NitscheStrategyPoro::complete_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, Teuchos::RCP<Core::LinAlg::SparseMatrix> kc)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_porofluid:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix())
              .GlobalAssemble(
                  *Global::Problem::instance()->get_dis("porofluid")->dof_row_map(),  // col map
                  *Global::Problem::instance()->get_dis("structure")->dof_row_map(),  // row map
                  true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::porofluid_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix())
              .GlobalAssemble(
                  *Global::Problem::instance()->get_dis("structure")->dof_row_map(),  // col map
                  *Global::Problem::instance()->get_dis("porofluid")->dof_row_map(),  // row map
                  true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::porofluid_porofluid:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix()).GlobalAssemble(true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::NitscheStrategy::complete_matrix_block_ptr(bt, kc);
      break;
  }
}

Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategyPoro::get_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bp) const
{
  if (!curr_state_eval_) FOUR_C_THROW("you didn't evaluate this contact state first");

  switch (bp)
  {
    case CONTACT::MatBlockType::porofluid_porofluid:
      return kpp_;
    case CONTACT::MatBlockType::porofluid_displ:
      return kpd_;
    case CONTACT::MatBlockType::displ_porofluid:
      return kdp_;
    default:
      return CONTACT::NitscheStrategy::get_matrix_block_ptr(bp, nullptr);
  }
}

FOUR_C_NAMESPACE_CLOSE
