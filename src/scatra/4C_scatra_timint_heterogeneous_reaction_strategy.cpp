/*----------------------------------------------------------------------*/
/*! \file

 \brief Solution strategy for heterogeneous reactions. This is not meshtying!!!

  \level 3

*/
/*----------------------------------------------------------------------*/
#include "4C_scatra_timint_heterogeneous_reaction_strategy.hpp"

#include "4C_discretization_dofset_gidbased_wrapper.hpp"
#include "4C_discretization_dofset_merged_wrapper.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_rebalance_print.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               vuong 06/16 |
 *----------------------------------------------------------------------*/
SCATRA::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy(
    SCATRA::ScaTraTimIntImpl* scatratimint)
    : MeshtyingStrategyStd(scatratimint), issetup_(false), isinit_(false)
{
  return;
}  // SCATRA::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy


/*------------------------------------------------------------------------*
 | evaluate heterogeneous reactions (actually no mesh tying    vuong 06/16 |
 *------------------------------------------------------------------------*/
void SCATRA::HeterogeneousReactionStrategy::EvaluateMeshtying()
{
  check_is_init();
  check_is_setup();

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_heteroreac_mat_and_rhs, condparams);

  // set global state vectors according to time-integration scheme
  discret_->set_state("phinp", scatratimint_->Phiafnp());
  discret_->set_state("hist", scatratimint_->Hist());

  // provide scatra discretization with convective velocity
  discret_->set_state(scatratimint_->NdsVel(), "convective velocity field",
      scatratimint_->discretization()->GetState(
          scatratimint_->NdsVel(), "convective velocity field"));

  // provide scatra discretization with velocity
  discret_->set_state(scatratimint_->NdsVel(), "velocity field",
      scatratimint_->discretization()->GetState(scatratimint_->NdsVel(), "velocity field"));

  if (scatratimint_->IsALE())
  {
    discret_->set_state(scatratimint_->NdsDisp(), "dispnp",
        scatratimint_->discretization()->GetState(scatratimint_->NdsDisp(), "dispnp"));
  }

  discret_->Evaluate(condparams, scatratimint_->SystemMatrix(), scatratimint_->Residual());

  // now we clear all states.
  // it would be nicer to do this directly before all
  // states are set at the beginning of this method.
  // However, in this case we are not able to set states externally
  // before this method is called.
  // See the call hierarchy of HeterogeneousReactionStrategy::set_state()
  // to check in which algorithms states are set on discret_ .
  discret_->ClearState();
  return;
}  // SCATRA::HeterogeneousReactionStrategy::EvaluateMeshtying


/*----------------------------------------------------------------------*
 | initialize meshtying objects                              rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::HeterogeneousReactionStrategy::setup_meshtying()
{
  // call Init() of base class
  SCATRA::MeshtyingStrategyStd::setup_meshtying();

  // make sure we set up everything properly
  heterogeneous_reaction_sanity_check();

  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(scatratimint_->discretization()->Comm().Clone());

  // standard case
  discret_ = Teuchos::rcp(new DRT::Discretization(
      scatratimint_->discretization()->Name(), com, GLOBAL::Problem::Instance()->NDim()));

  // call complete without assigning degrees of freedom
  discret_->fill_complete(false, true, false);

  Teuchos::RCP<DRT::Discretization> scatradis = scatratimint_->discretization();

  // create scatra elements if the scatra discretization is empty
  {
    // fill scatra discretization by cloning fluid discretization
    CORE::FE::CloneDiscretizationFromCondition<SCATRA::ScatraReactionCloneStrategy>(*scatradis,
        *discret_, "ScatraHeteroReactionSlave", GLOBAL::Problem::Instance()->CloningMaterialMap());

    // set implementation type of cloned scatra elements
    for (int i = 0; i < discret_->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::Transport* element =
          dynamic_cast<DRT::ELEMENTS::Transport*>(discret_->lColElement(i));
      if (element == nullptr) FOUR_C_THROW("Invalid element type!");

      if (element->Material()->MaterialType() == CORE::Materials::m_matlist_reactions)
        element->SetImplType(INPAR::SCATRA::impltype_advreac);
      else
        FOUR_C_THROW("Invalid material type for HeterogeneousReactionStrategy!");
    }  // loop over all column elements
  }

  {
    // build a dofset that merges the DOFs from both sides
    // convention: the order of the merged dofset will be
    //  _            _
    // | slave dofs   |
    // |              |
    // |_master dofs _|
    //
    // slave side is supposed to be the surface discretization
    //
    Teuchos::RCP<CORE::Dofsets::DofSetMergedWrapper> newdofset =
        Teuchos::rcp(new CORE::Dofsets::DofSetMergedWrapper(scatradis->GetDofSetProxy(), scatradis,
            "ScatraHeteroReactionMaster", "ScatraHeteroReactionSlave"));

    // assign the dofset to the reaction discretization
    discret_->ReplaceDofSet(newdofset, false);

    // add all secondary dofsets as proxies
    for (int ndofset = 1; ndofset < scatratimint_->discretization()->NumDofSets(); ++ndofset)
    {
      Teuchos::RCP<CORE::Dofsets::DofSetGIDBasedWrapper> gidmatchingdofset =
          Teuchos::rcp(new CORE::Dofsets::DofSetGIDBasedWrapper(scatratimint_->discretization(),
              scatratimint_->discretization()->GetDofSetProxy(ndofset)));
      discret_->AddDofSet(gidmatchingdofset);
    }

    // done. Rebuild all maps and boundary condition geometries
    discret_->fill_complete(true, true, true);

    if (com->MyPID() == 0 and com->NumProc() > 1)
      std::cout << "parallel distribution of auxiliary discr. with standard ghosting" << std::endl;
    CORE::REBALANCE::UTILS::print_parallel_distribution(*discret_);
  }

  set_is_setup(true);
  return;
}


/*----------------------------------------------------------------------*
 | setup meshtying objects                                  vuong 06/16 |
 *----------------------------------------------------------------------*/
void SCATRA::HeterogeneousReactionStrategy::InitMeshtying()
{
  set_is_setup(false);

  // call Init() of base class
  SCATRA::MeshtyingStrategyStd::InitMeshtying();

  set_is_init(true);
  return;
}


/*----------------------------------------------------------------------*
 | Evaluate conditioned elements                            rauch 08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::HeterogeneousReactionStrategy::evaluate_condition(Teuchos::ParameterList& params,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3, const std::string& condstring, const int condid)
{
  check_is_init();
  check_is_setup();

  // Call evaluate_condition on auxiliary discretization.
  // This condition has all dofs, both from the volume-
  // bound scalars and from the surface-bound scalars.
  discret_->evaluate_condition(params, systemmatrix1, systemmatrix2, systemvector1, systemvector2,
      systemvector3, condstring, condid);

  return;
}


/*----------------------------------------------------------------------*
 | Set state on auxiliary discretization                    rauch 12/16 |
 *----------------------------------------------------------------------*/
void SCATRA::HeterogeneousReactionStrategy::set_state(
    unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state)
{
  discret_->set_state(nds, name, state);
  return;
}


/*----------------------------------------------------------------------*
 | sanity check for some assumptions and conventions        rauch 06/17 |
 *----------------------------------------------------------------------*/
void SCATRA::HeterogeneousReactionStrategy::heterogeneous_reaction_sanity_check()
{
  bool valid_slave = false;

  const Epetra_Comm& com = scatratimint_->discretization()->Comm();

  if (com.MyPID() == 0) std::cout << " Sanity check for HeterogeneousReactionStrategy ...";

  CORE::Conditions::Condition* slave_cond =
      scatratimint_->discretization()->GetCondition("ScatraHeteroReactionSlave");

  const Epetra_Map* element_row_map = scatratimint_->discretization()->ElementRowMap();

  // loop over row elements
  for (int lid = 0; lid < element_row_map->NumMyElements(); lid++)
  {
    const int gid = element_row_map->GID(lid);

    CORE::Elements::Element* ele = scatratimint_->discretization()->gElement(gid);
    CORE::Nodes::Node** nodes = ele->Nodes();
    if (ele->Shape() == CORE::FE::CellType::quad4 or ele->Shape() == CORE::FE::CellType::tri3)
    {
      for (int node = 0; node < ele->num_node(); node++)
      {
        const int node_gid = nodes[node]->Id();

        if (not slave_cond->ContainsNode(node_gid))
        {
          FOUR_C_THROW(
              "Surface discretization for membrane transport is "
              "supposed to wear ScatraHeteroReactionSlave condition!");
        }
        else
        {
          valid_slave = true;
          break;
        }
      }  // loop over nodes of row ele

      if (valid_slave) break;
    }  // if surface transport element

    else if (ele->Shape() == CORE::FE::CellType::hex8 or ele->Shape() == CORE::FE::CellType::tet4)
    {
      // no check so far
    }  // if volume transport element

    else
    {
      FOUR_C_THROW(
          "please implement check for new combination of volume transport "
          "- surface transport elements.");
    }

  }  // loop over row elements


  com.Barrier();
  if (com.MyPID() == 0) std::cout << " Passed." << std::endl;

  return;
}

FOUR_C_NAMESPACE_CLOSE
