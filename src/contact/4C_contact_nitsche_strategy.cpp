/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_nitsche_strategy.hpp"

#include "4C_contact_interface.hpp"
#include "4C_contact_nitsche_utils.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_so3_plast_ssn.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | global evaluation method called from time integrator     seitz 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::NitscheStrategy::apply_force_stiff_cmt(Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
    const int iter, bool predictor)
{
  // mortar initialization and evaluation
  set_state(Mortar::state_new_displacement, *dis);

  // just a Nitsche-version
  Teuchos::RCP<Epetra_FEVector> fc = Teuchos::rcp(new Epetra_FEVector(f->Map()));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> kc = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      (dynamic_cast<Epetra_CrsMatrix*>(&(*kt->epetra_operator())))->RowMap(), 100, true, false,
      Core::LinAlg::SparseMatrix::FE_MATRIX));

  // Evaluation for all interfaces
  for (const auto& interface : interface_)
  {
    interface->initialize();
    interface->evaluate(0, step_, iter_);
    for (int e = 0; e < interface->discret().element_col_map()->NumMyElements(); ++e)
    {
      auto* mele = dynamic_cast<Mortar::Element*>(
          interface->discret().g_element(interface->discret().element_col_map()->GID(e)));
      mele->get_nitsche_container().assemble_rhs(mele, CONTACT::VecBlockType::displ, fc);
      mele->get_nitsche_container().assemble_matrix(mele, CONTACT::MatBlockType::displ_displ, kc);
    }
  }

  // now we also did this state
  curr_state_eval_ = true;

  if (fc->GlobalAssemble(Add, false) != 0) FOUR_C_THROW("GlobalAssemble failed");
  // add negative contact force here since the time integrator handed me a rhs!
  if (f->Update(-1., *fc, 1.)) FOUR_C_THROW("update went wrong");
  dynamic_cast<Epetra_FECrsMatrix&>(*kc->epetra_matrix()).GlobalAssemble(true, Add);
  kt->un_complete();
  kt->add(*kc, false, 1., 1.);
  kt->complete();
}


/*----------------------------------------------------------------------*
 |  read restart information for contact                     seitz 10/16|
 *----------------------------------------------------------------------*/
void CONTACT::NitscheStrategy::do_read_restart(Core::IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // check whether this is a restart with contact of a previously
  // non-contact simulation run (if yes, we have to be careful not
  // to try to read certain, in this case non-existing, vectors
  // such as the activetoggle or sliptoggle vectors, but rather
  // initialize the restart active and slip sets as being empty)
  bool restartwithcontact = Core::UTILS::IntegralValue<int>(params(), "RESTART_WITH_CONTACT");
  if (restartwithcontact) FOUR_C_THROW("not supported for nitsche contact");

  // set restart displacement state
  set_state(Mortar::state_new_displacement, *dis);
  set_state(Mortar::state_old_displacement, *dis);

  // Evaluation for all interfaces
  for (const auto& interface : interface_) interface->initialize();

  if (friction_)
  {
    for (const auto& interface : interface_)
    {
      interface->evaluate_nodal_normals();
      interface->export_nodal_normals();
    }
    store_to_old(Mortar::StrategyBase::n_old);
  }

  if (Core::UTILS::IntegralValue<int>(params(), "NITSCHE_PENALTY_ADAPTIVE"))
    update_trace_ineq_etimates();
}

void CONTACT::NitscheStrategy::set_state(
    const enum Mortar::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == Mortar::state_new_displacement)
  {
    double inf_delta = 0.;
    if (curr_state_ == Teuchos::null)
    {
      curr_state_ = Teuchos::rcp(new Epetra_Vector(vec));
      inf_delta = 1.e12;
    }
    else
    {
      Epetra_Vector delta(vec);
      delta.Update(-1., *curr_state_, 1.);
      delta.NormInf(&inf_delta);
    }
    if (inf_delta < 1.e-12)
      return;
    else
    {
      curr_state_eval_ = false;
      (*curr_state_) = vec;
      AbstractStrategy::set_state(statename, vec);
      set_parent_state(statename, vec);
    }
  }
  else
  {
    curr_state_eval_ = false;
    AbstractStrategy::set_state(statename, vec);
  }
}

/*------------------------------------------------------------------------*
 |                                                             seitz 10/16|
 *------------------------------------------------------------------------*/
void CONTACT::NitscheStrategy::set_parent_state(
    const enum Mortar::StateType& statename, const Epetra_Vector& vec)
{
  Teuchos::RCP<Core::FE::Discretization> dis = Global::Problem::instance()->get_dis("structure");
  if (dis == Teuchos::null) FOUR_C_THROW("didn't get my discretization");
  if (statename == Mortar::state_new_displacement || statename == Mortar::state_svelocity)
  {
    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*dis->dof_col_map(), true));
    Core::LinAlg::Export(vec, *global);

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

        // this gets values in local order
        ele->parent_element()->location_vector(*dis, lm, lmowner, lmstride);

        std::vector<double> myval;
        Core::FE::ExtractMyValues(*global, myval, lm);

        switch (statename)
        {
          case Mortar::state_new_displacement:
          {
            ele->mo_data().parent_disp() = myval;
            ele->mo_data().parent_dof() = lm;
            break;
          }
          case Mortar::state_svelocity:
          {
            ele->mo_data().parent_vel() = myval;
            break;
          }
          default:
            FOUR_C_THROW("unknown statename");
        }
      }
    }
  }
}

void CONTACT::NitscheStrategy::evaluate_force(CONTACT::ParamsInterface& cparams)
{
  integrate(cparams);
}

void CONTACT::NitscheStrategy::evaluate_force_stiff(CONTACT::ParamsInterface& cparams)
{
  integrate(cparams);
}

void CONTACT::NitscheStrategy::reset(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp, const Epetra_Vector& xnew)
{
  set_state(Mortar::state_new_displacement, dispnp);
}

void CONTACT::NitscheStrategy::run_post_compute_x(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  // do nothing
}

void CONTACT::NitscheStrategy::integrate(const CONTACT::ParamsInterface& cparams)
{
  // we already did this displacement state
  if (curr_state_eval_) return;

  // time measurement (on each processor)
  const double t_start = Teuchos::Time::wallTime();

  // Evaluation for all interfaces
  for (const auto& interface : interface_)
  {
    interface->interface_params().set<double>("TIMESTEP", cparams.get_delta_time());
    interface->initialize();
    interface->evaluate(0, step_, iter_);

    // store required integration time
    inttime_ += interface->inttime();
  }

  // check the parallel distribution
  check_parallel_distribution(t_start);

  // now we also did this state
  curr_state_eval_ = true;

  // ... and we can assemble the matrix and rhs
  fc_ = create_rhs_block_ptr(CONTACT::VecBlockType::displ);
  kc_ = create_matrix_block_ptr(CONTACT::MatBlockType::displ_displ);
}

Teuchos::RCP<Epetra_FEVector> CONTACT::NitscheStrategy::setup_rhs_block_vec(
    const enum CONTACT::VecBlockType& bt) const
{
  switch (bt)
  {
    case CONTACT::VecBlockType::displ:
      return Teuchos::rcp(
          new Epetra_FEVector(*Global::Problem::instance()->get_dis("structure")->dof_row_map()));
    default:
      FOUR_C_THROW("you should not be here");
      break;
  }
  return Teuchos::null;
}

Teuchos::RCP<Epetra_FEVector> CONTACT::NitscheStrategy::create_rhs_block_ptr(
    const enum CONTACT::VecBlockType& bt) const
{
  if (!curr_state_eval_) FOUR_C_THROW("you didn't evaluate this contact state first");

  Teuchos::RCP<Epetra_FEVector> fc = setup_rhs_block_vec(bt);

  for (const auto& interface : interface_)
  {
    for (int e = 0; e < interface->discret().element_col_map()->NumMyElements(); ++e)
    {
      auto* mele = dynamic_cast<Mortar::Element*>(
          interface->discret().g_element(interface->discret().element_col_map()->GID(e)));
      auto& nitsche_container = mele->get_nitsche_container();
      nitsche_container.assemble_rhs(mele, bt, fc);
    }
  }
  if (fc->GlobalAssemble(Add, false) != 0) FOUR_C_THROW("GlobalAssemble failed");

  return fc;
}

Teuchos::RCP<const Epetra_Vector> CONTACT::NitscheStrategy::get_rhs_block_ptr(
    const enum CONTACT::VecBlockType& bt) const
{
  if (!curr_state_eval_)
    FOUR_C_THROW("you didn't evaluate this contact state for %s first",
        CONTACT::VecBlockTypeToStr(bt).c_str());

  switch (bt)
  {
    case CONTACT::VecBlockType::displ:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(fc_), 0));
    case CONTACT::VecBlockType::constraint:
      return Teuchos::null;
    default:
      FOUR_C_THROW("get_rhs_block_ptr: your type is no treated properly!");
      break;
  }

  return Teuchos::null;
}

Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategy::setup_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_displ:
      return Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *Global::Problem::instance()->get_dis("structure")->dof_row_map()),
          100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    default:
      FOUR_C_THROW("you should not be here");
      break;
  }
  return Teuchos::null;
}

void CONTACT::NitscheStrategy::complete_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, Teuchos::RCP<Core::LinAlg::SparseMatrix> kc)
{
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_displ:
      kc->complete();
      break;
    default:
      FOUR_C_THROW("you should not be here");
      break;
  }
}

Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategy::create_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt)
{
  if (!curr_state_eval_) FOUR_C_THROW("you didn't evaluate this contact state first");

  Teuchos::RCP<Core::LinAlg::SparseMatrix> kc = setup_matrix_block_ptr(bt);

  for (const auto& interface : interface_)
  {
    for (int e = 0; e < interface->discret().element_col_map()->NumMyElements(); ++e)
    {
      auto* mele = dynamic_cast<Mortar::Element*>(
          interface->discret().g_element(interface->discret().element_col_map()->GID(e)));
      mele->get_nitsche_container().assemble_matrix(mele, bt, kc);
    }
  }

  complete_matrix_block_ptr(bt, kc);

  return kc;
}

Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::NitscheStrategy::get_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  if (!curr_state_eval_) FOUR_C_THROW("you didn't evaluate this contact state first");

  if (bt == CONTACT::MatBlockType::displ_displ)
    return kc_;
  else
    FOUR_C_THROW("get_matrix_block_ptr: your type is no treated properly!");

  return Teuchos::null;
}

void CONTACT::NitscheStrategy::setup(bool redistributed, bool init)
{
  // we need to init the isselfcontact_ flag here, as we do not want to call the AbstractStrategy
  if (init)
  {
    // set potential global self contact status
    // (this is TRUE if at least one contact interface is a self contact interface)
    bool selfcontact = false;
    for (const auto& interface : interfaces())
      if (interface->self_contact()) selfcontact = true;

    if (selfcontact) isselfcontact_ = true;
  }
  reconnect_parent_elements();
  curr_state_ = Teuchos::null;
  curr_state_eval_ = false;
}

void CONTACT::NitscheStrategy::update_trace_ineq_etimates()
{
  auto NitWgt =
      Core::UTILS::IntegralValue<Inpar::CONTACT::NitscheWeighting>(params(), "NITSCHE_WEIGHTING");
  for (const auto& interface : interface_)
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

void CONTACT::NitscheStrategy::update(Teuchos::RCP<const Epetra_Vector> dis)
{
  if (Core::UTILS::IntegralValue<int>(params(), "NITSCHE_PENALTY_ADAPTIVE"))
    update_trace_ineq_etimates();
  if (friction_)
  {
    store_to_old(Mortar::StrategyBase::n_old);
    set_state(Mortar::state_old_displacement, *dis);
  }
}

void CONTACT::NitscheStrategy::evaluate_reference_state()
{
  if (friction_)
  {
    for (const auto& interface : interface_)
    {
      interface->evaluate_nodal_normals();
      interface->export_nodal_normals();
    }
    store_to_old(Mortar::StrategyBase::n_old);
  }

  update_trace_ineq_etimates();
}


/*----------------------------------------------------------------------------------------------*
 |  Reconnect Contact Element -- Parent Element Pointers (required for restart)       ager 04/16|
 *---------------------------------------------------------------------------------------------*/
void CONTACT::NitscheStrategy::reconnect_parent_elements()
{
  Teuchos::RCP<Core::FE::Discretization> voldis = Global::Problem::instance()->get_dis("structure");

  for (const auto& contact_interface : contact_interfaces())
  {
    const Epetra_Map* elecolmap = voldis->element_col_map();

    const Epetra_Map* ielecolmap = contact_interface->discret().element_col_map();

    for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
    {
      const int gid = ielecolmap->GID(i);

      Core::Elements::Element* ele = contact_interface->discret().g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
      auto* faceele = dynamic_cast<Core::Elements::FaceElement*>(ele);

      const int volgid = faceele->parent_element_id();
      if (elecolmap->LID(volgid) == -1)  // Volume discretization has not Element
        FOUR_C_THROW(
            "Manager::reconnect_parent_elements: Element %d does not exist on this Proc!", volgid);

      Core::Elements::Element* vele = voldis->g_element(volgid);
      if (!vele) FOUR_C_THROW("Cannot find element with gid %", volgid);

      faceele->set_parent_master_element(vele, faceele->face_parent_number());

      auto* vele_plast = dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>*>(vele);
      if (vele_plast) vele_plast->set_is_nitsche_contact_ele(true);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
