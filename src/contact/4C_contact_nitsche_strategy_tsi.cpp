/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_nitsche_strategy_tsi.hpp"

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

void CONTACT::NitscheStrategyTsi::set_state(
    const enum Mortar::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == Mortar::state_temperature)
  {
    double inf_delta = 0.;
    if (curr_state_temp_ == Teuchos::null)
    {
      curr_state_temp_ = Teuchos::rcp(new Epetra_Vector(vec));
      inf_delta = 1.e12;
    }
    else
    {
      Epetra_Vector delta(vec);
      delta.Update(-1., *curr_state_temp_, 1.);
      delta.NormInf(&inf_delta);
    }
    if (inf_delta < 1.e-16)
      return;
    else
    {
      curr_state_eval_ = false;
      (*curr_state_temp_) = vec;
      set_parent_state(statename, vec);
    }
  }
  else
    CONTACT::NitscheStrategy::set_state(statename, vec);
}

/*------------------------------------------------------------------------*
 |                                                             seitz 10/16|
 *------------------------------------------------------------------------*/
void CONTACT::NitscheStrategyTsi::set_parent_state(
    const enum Mortar::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == Mortar::state_temperature)
  {
    Teuchos::RCP<Core::FE::Discretization> disT = Global::Problem::instance()->get_dis("thermo");
    if (disT.is_null()) FOUR_C_THROW("set state temperature, but no thermo-discretization???");

    Teuchos::RCP<Epetra_Vector> global =
        Teuchos::rcp(new Epetra_Vector(*disT->dof_col_map(), true));
    Core::LinAlg::Export(vec, *global);

    // set state on interfaces
    for (auto& interface : interface_)
    {
      Core::FE::Discretization& idiscret = interface->discret();

      for (int j = 0; j < idiscret.element_col_map()->NumMyElements(); ++j)
      {
        int gid = idiscret.element_col_map()->GID(j);

        Core::Elements::Element* e = idiscret.g_element(gid);
        if (e == nullptr) FOUR_C_THROW("basic element not found");

        auto* ele = dynamic_cast<Mortar::Element*>(idiscret.g_element(gid));
        Core::Elements::Element* ele_parentT = disT->g_element(ele->parent_element_id());

        std::vector<int> lm, lmowner, lmstride;
        ele_parentT->location_vector(*disT, lm, lmowner, lmstride);

        std::vector<double> myval;
        Core::FE::ExtractMyValues(*global, myval, lm);

        ele->mo_data().parent_temp() = myval;
        ele->mo_data().parent_temp_dof() = lm;
      }
    }
  }
  else
    CONTACT::NitscheStrategy::set_parent_state(statename, vec);
}

void CONTACT::NitscheStrategyTsi::setup(bool redistributed, bool init)
{
  CONTACT::NitscheStrategy::setup(redistributed, init);

  curr_state_temp_ = Teuchos::null;
}

void CONTACT::NitscheStrategyTsi::update_trace_ineq_etimates()
{
  auto NitWgt =
      Core::UTILS::IntegralValue<Inpar::CONTACT::NitscheWeighting>(params(), "NITSCHE_WEIGHTING");
  for (auto& interface : interface_)
  {
    for (int e = 0; e < interface->discret().element_col_map()->NumMyElements(); ++e)
    {
      auto* mele = dynamic_cast<Mortar::Element*>(
          interface->discret().g_element(interface->discret().element_col_map()->GID(e)));
      if (NitWgt == Inpar::CONTACT::NitWgt_slave && !mele->is_slave()) continue;
      if (NitWgt == Inpar::CONTACT::NitWgt_master && mele->is_slave()) continue;
      mele->estimate_nitsche_trace_max_eigenvalue_combined();
    }
  }
}

Teuchos::RCP<Epetra_FEVector> CONTACT::NitscheStrategyTsi::setup_rhs_block_vec(
    const enum CONTACT::VecBlockType& bt) const
{
  switch (bt)
  {
    case CONTACT::VecBlockType::temp:
      return Teuchos::rcp(
          new Epetra_FEVector(*Global::Problem::instance()->get_dis("thermo")->dof_row_map()));
    default:
      return CONTACT::NitscheStrategy::setup_rhs_block_vec(bt);
  }
}

Teuchos::RCP<const Epetra_Vector> CONTACT::NitscheStrategyTsi::get_rhs_block_ptr(
    const enum CONTACT::VecBlockType& bt) const
{
  if (bt == CONTACT::VecBlockType::constraint) return Teuchos::null;

  if (!curr_state_eval_)
    FOUR_C_THROW(
        "you didn't evaluate this contact state for %s first", VecBlockTypeToStr(bt).c_str());

  switch (bt)
  {
    case CONTACT::VecBlockType::temp:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(ft_), 0));
    default:
      return CONTACT::NitscheStrategy::get_rhs_block_ptr(bt);
  }
}

Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategyTsi::setup_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_temp:
      return Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *Global::Problem::instance()->get_dis("structure")->dof_row_map()),
          100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    case CONTACT::MatBlockType::temp_displ:
    case CONTACT::MatBlockType::temp_temp:
      return Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *Global::Problem::instance()->get_dis("thermo")->dof_row_map()),
          100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    default:
      return CONTACT::NitscheStrategy::setup_matrix_block_ptr(bt);
  }
}

void CONTACT::NitscheStrategyTsi::complete_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, Teuchos::RCP<Core::LinAlg::SparseMatrix> kc)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_temp:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix())
              .GlobalAssemble(
                  *Global::Problem::instance()->get_dis("thermo")->dof_row_map(),     // col map
                  *Global::Problem::instance()->get_dis("structure")->dof_row_map(),  // row map
                  true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::temp_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix())
              .GlobalAssemble(
                  *Global::Problem::instance()->get_dis("structure")->dof_row_map(),  // col map
                  *Global::Problem::instance()->get_dis("thermo")->dof_row_map(),     // row map
                  true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::temp_temp:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix()).GlobalAssemble(true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::NitscheStrategy::complete_matrix_block_ptr(bt, kc);
      break;
  }
}

Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategyTsi::get_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, const ParamsInterface* cparams) const
{
  if (!curr_state_eval_) FOUR_C_THROW("you didn't evaluate this contact state first");

  switch (bt)
  {
    case CONTACT::MatBlockType::temp_temp:
      return ktt_;
    case CONTACT::MatBlockType::temp_displ:
      return ktd_;
    case CONTACT::MatBlockType::displ_temp:
      return kdt_;
    default:
      return CONTACT::NitscheStrategy::get_matrix_block_ptr(bt, cparams);
  }
}


void CONTACT::NitscheStrategyTsi::integrate(const CONTACT::ParamsInterface& cparams)
{
  CONTACT::NitscheStrategy::integrate(cparams);

  ft_ = create_rhs_block_ptr(CONTACT::VecBlockType::temp);
  ktt_ = create_matrix_block_ptr(CONTACT::MatBlockType::temp_temp);
  ktd_ = create_matrix_block_ptr(CONTACT::MatBlockType::temp_displ);
  kdt_ = create_matrix_block_ptr(CONTACT::MatBlockType::displ_temp);
}

FOUR_C_NAMESPACE_CLOSE
