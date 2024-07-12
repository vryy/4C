/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "4C_contact_nitsche_strategy_ssi.hpp"

#include "4C_contact_interface.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_element.hpp"

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::NitscheStrategySsi::integrate(const CONTACT::ParamsInterface& cparams)
{
  CONTACT::NitscheStrategy::integrate(cparams);

  fs_ = create_rhs_block_ptr(CONTACT::VecBlockType::scatra);
  kss_ = create_matrix_block_ptr(CONTACT::MatBlockType::scatra_scatra);
  ksd_ = create_matrix_block_ptr(CONTACT::MatBlockType::scatra_displ);
  kds_ = create_matrix_block_ptr(CONTACT::MatBlockType::displ_scatra);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::NitscheStrategySsi::evaluate_reference_state()
{
  // initialize an estimate of TraceHE
  init_trace_he();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::NitscheStrategySsi::init_trace_he()
{
  for (const auto& interface : interface_)
  {
    for (int e = 0; e < interface->discret().element_col_map()->NumMyElements(); ++e)
    {
      auto* mele = dynamic_cast<Mortar::Element*>(
          interface->discret().g_element(interface->discret().element_col_map()->GID(e)));
      // set the initial TraceHE to the edge length size of the element
      mele->trace_he() = mele->max_edge_size();
    }
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::NitscheStrategySsi::set_state(
    const enum Mortar::StateType& statename, const Epetra_Vector& vec)
{
  switch (statename)
  {
    case Mortar::state_elch:
    case Mortar::state_scalar:
    {
      double inf_delta = 0.0;
      if (curr_state_scalar_ == Teuchos::null)
      {
        curr_state_scalar_ = Teuchos::rcp(new Epetra_Vector(vec));
        inf_delta = 1.0e12;
      }
      else
      {
        auto delta(vec);
        delta.Update(-1.0, *curr_state_scalar_, 1.0);
        delta.NormInf(&inf_delta);
      }
      if (inf_delta < 1.0e-16)
        return;
      else
      {
        curr_state_eval_ = false;
        *curr_state_scalar_ = vec;
        set_parent_state(statename, vec);
      }
      break;
    }
    default:
    {
      CONTACT::NitscheStrategy::set_state(statename, vec);
      break;
    }
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::NitscheStrategySsi::set_parent_state(
    const enum Mortar::StateType& statename, const Epetra_Vector& vec)
{
  switch (statename)
  {
    case Mortar::state_elch:
    case Mortar::state_scalar:
    {
      auto scatra_dis = Global::Problem::instance()->get_dis("scatra");
      if (scatra_dis == Teuchos::null) FOUR_C_THROW("didn't get scatra discretization");

      auto scatra_dofcolmap = Teuchos::rcp(new Epetra_Vector(*scatra_dis->dof_col_map(), true));
      Core::LinAlg::export_to(vec, *scatra_dofcolmap);

      // set state on interfaces
      for (const auto& interface : interface_)
      {
        // get the interface discretization
        const auto& interface_dis = interface->discret();

        // loop over all interface column elements owned by current processor
        for (int j = 0; j < interface_dis.element_col_map()->NumMyElements(); ++j)
        {
          const int gid = interface_dis.element_col_map()->GID(j);

          auto* interface_ele = interface_dis.g_element(gid);
          if (interface_ele == nullptr) FOUR_C_THROW("Did not find element.");

          auto* mortar_ele = dynamic_cast<Mortar::Element*>(interface_ele);
          auto* mortar_parent_ele = scatra_dis->g_element(mortar_ele->parent_element_id());

          std::vector<int> lm;
          std::vector<int> lmowner;
          std::vector<int> lmstride;

          if (mortar_parent_ele == nullptr)
            FOUR_C_THROW("Did not get parent element to extract scalar values");
          else
            mortar_parent_ele->location_vector(*scatra_dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          Core::FE::ExtractMyValues(*scatra_dofcolmap, myval, lm);

          mortar_ele->mo_data().parent_scalar() = myval;
          mortar_ele->mo_data().parent_scalar_dof() = lm;
        }
      }
      break;
    }
    default:
    {
      CONTACT::NitscheStrategy::set_parent_state(statename, vec);
      break;
    }
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_FEVector> CONTACT::NitscheStrategySsi::setup_rhs_block_vec(
    const enum CONTACT::VecBlockType& bt) const
{
  switch (bt)
  {
    case CONTACT::VecBlockType::elch:
    case CONTACT::VecBlockType::scatra:
      return Teuchos::rcp(
          new Epetra_FEVector(*Global::Problem::instance()->get_dis("scatra")->dof_row_map()));
    default:
      return CONTACT::NitscheStrategy::setup_rhs_block_vec(bt);
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::NitscheStrategySsi::get_rhs_block_ptr(
    const enum CONTACT::VecBlockType& bp) const
{
  if (bp == CONTACT::VecBlockType::constraint) return Teuchos::null;

  if (!curr_state_eval_)
    FOUR_C_THROW("you didn't evaluate this contact state for %s first",
        CONTACT::VecBlockTypeToStr(bp).c_str());

  switch (bp)
  {
    case CONTACT::VecBlockType::elch:
    case CONTACT::VecBlockType::scatra:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(fs_), 0));
    default:
      return CONTACT::NitscheStrategy::get_rhs_block_ptr(bp);
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategySsi::setup_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_elch:
    case CONTACT::MatBlockType::displ_scatra:
      return Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *Global::Problem::instance()->get_dis("structure")->dof_row_map()),
          100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    case CONTACT::MatBlockType::elch_displ:
    case CONTACT::MatBlockType::elch_elch:
    case CONTACT::MatBlockType::scatra_displ:
    case CONTACT::MatBlockType::scatra_scatra:
      return Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *Global::Problem::instance()->get_dis("scatra")->dof_row_map()),
          100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    default:
      return CONTACT::NitscheStrategy::setup_matrix_block_ptr(bt);
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::NitscheStrategySsi::complete_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, Teuchos::RCP<Core::LinAlg::SparseMatrix> kc)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_elch:
    case CONTACT::MatBlockType::displ_scatra:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix())
              .GlobalAssemble(
                  *Global::Problem::instance()->get_dis("scatra")->dof_row_map(),     // col map
                  *Global::Problem::instance()->get_dis("structure")->dof_row_map(),  // row map
                  true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::elch_displ:
    case CONTACT::MatBlockType::scatra_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix())
              .GlobalAssemble(
                  *Global::Problem::instance()->get_dis("structure")->dof_row_map(),  // col map
                  *Global::Problem::instance()->get_dis("scatra")->dof_row_map(),     // row map
                  true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    case CONTACT::MatBlockType::elch_elch:
    case CONTACT::MatBlockType::scatra_scatra:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix()).GlobalAssemble(true, Add))
        FOUR_C_THROW("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::NitscheStrategy::complete_matrix_block_ptr(bt, kc);
      break;
  }
}

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategySsi::get_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bp) const
{
  if (!curr_state_eval_) FOUR_C_THROW("you didn't evaluate this contact state first");

  switch (bp)
  {
    case CONTACT::MatBlockType::elch_elch:
    case CONTACT::MatBlockType::scatra_scatra:
      return kss_;
    case CONTACT::MatBlockType::elch_displ:
    case CONTACT::MatBlockType::scatra_displ:
      return ksd_;
    case CONTACT::MatBlockType::displ_elch:
    case CONTACT::MatBlockType::displ_scatra:
      return kds_;
    default:
      return CONTACT::NitscheStrategy::get_matrix_block_ptr(bp, nullptr);
  }
}

FOUR_C_NAMESPACE_CLOSE
