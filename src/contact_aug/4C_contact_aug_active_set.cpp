/*----------------------------------------------------------------------------*/
/*! \file
\brief Identify the correct active set.


\level 2
*/
/*----------------------------------------------------------------------------*/

#include "4C_contact_aug_active_set.hpp"

#include "4C_contact_aug_interface.hpp"
#include "4C_contact_aug_potential.hpp"
#include "4C_contact_aug_strategy.hpp"
#include "4C_contact_node.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::ActiveSet::Update(const CONTACT::ParamsInterface& cparams)
{
  if (skip_update()) return;

  const Status gstatus = update_status(cparams);

  update_maps(cparams);

  post_update(cparams, gstatus);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::Aug::ActiveSet::skip_update() const
{
  const DataContainer& data = strategy_.data();
  const plain_interface_set& interfaces = strategy_.interfaces();

  // get out of here if not in the semi-smooth Newton case
  // (but before doing this, check if there are invalid active nodes)
  if (not data.IsSemiSmoothNewton())
  {
    // loop over all interfaces
    for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end();
         ++cit)
    {
      const CONTACT::Interface& interface = **cit;

      // loop over all slave nodes on the current interface
      for (int j = 0; j < interface.SlaveRowNodes()->NumMyElements(); ++j)
      {
        int gid = interface.SlaveRowNodes()->GID(j);
        Core::Nodes::Node* node = interface.Discret().gNode(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* cnode = static_cast<Node*>(node);

        /* The nested active set strategy cannot deal with the case of
         * active nodes that have no integration segments/cells attached,
         * as this leads to zero rows in D and M and thus to singular systems.
         * However, this case might possibly happen when slave nodes slide
         * over the edge of a master body within one fixed active set step.
         * (Remark: Semi-smooth Newton has no problems in this case, as it
         * updates the active set after EACH Newton step, see below, and thus
         * would always set the corresponding nodes to INACTIVE.) */
        if (cnode->Active() && !cnode->HasSegment())
          FOUR_C_THROW("Active node %i without any segment/cell attached", cnode->Id());
      }
    }
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::ActiveSet::post_update(
    const CONTACT::ParamsInterface& cparams, const enum Status gstatus)
{
  DataContainer& data = strategy_.data();

  // check the convergence of the active set
  data.is_active_set_converged() = data.GActiveNodeRowMap().SameAs(data.g_old_active_slave_nodes());

  sanity_check(cparams, gstatus);

  if (gstatus == Status::changed)
  {
    data.SetVectorMapsValid(false);
    data.SetMatrixMapsValid(false);

    // reset the active/inactive state vectors
    data.Potential().Setup();
  }

  // set the new active/inactive state
  data.Potential().set_active_inactive_state();

  // update the history information only if it's no correction step of the active set
  if (cparams.IsDefaultStep())
  {
    // update flag for the contact status of the last iterate (history information)
    if (strategy_.IsInContact())
      data.was_in_contact_last_iter() = true;
    else
      data.was_in_contact_last_iter() = false;
  }
  // update flag for global contact status
  if (data.GActiveNodeRowMap().NumGlobalElements())
  {
    data.IsInContact() = true;
    data.WasInContact() = true;
  }
  else
    data.IsInContact() = false;

  //  int icount = 0;
  //  for ( plain_interface_set::const_iterator cit = interface_.begin();
  //        cit != interface_.end(); ++cit, ++icount )
  //  {
  //    Interface& interface = **cit;
  //    interface.write_nodal_coordinates_to_file( icount, *Data().g_active_node_row_map_ptr(),
  //        "../o/half_sphere/aug_nurbs_complete_active_slave_node_coordinates.data");
  //  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::ActiveSet::Status CONTACT::Aug::ActiveSet::update_status(
    const CONTACT::ParamsInterface& cparams) const
{
  plain_interface_set& interfaces = strategy_.interfaces();
  CONTACT::Aug::DataContainer& data = strategy_.data();

  // assume that active set has converged and check for opposite
  strategy_.data().is_active_set_converged() = true;

  std::vector<Status> istatus(interfaces.size(), Status::unchanged);

  // loop over all interfaces
  unsigned ilid = 0;
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end();
       ++cit, ++ilid)
  {
    const Interface& interface = dynamic_cast<const Interface&>(**cit);

    // loop over all slave nodes of the current interface
    const int num_my_slave_row_nodes = interface.SlaveRowNodes()->NumMyElements();
    int* my_slave_row_node_gids = interface.SlaveRowNodes()->MyGlobalElements();
    for (int j = 0; j < num_my_slave_row_nodes; ++j)
    {
      const int gid = my_slave_row_node_gids[j];
      Node* cnode = dynamic_cast<Node*>(interface.Discret().gNode(gid));
      if (!cnode) FOUR_C_THROW("Cannot find node with gid %", gid);

      /* read weighting factor cn
       * (this is necessary in semi-smooth Newton case, as the search for the
       * active set is now part of the Newton iteration. Thus, we do not know
       * the active / inactive status in advance and we can have a state in
       * which both the condition znormal = 0 and wgap = 0 are violated. Here
       * we have to weight the two violations via cn! */
      const int cn_lid = data.Cn().Map().LID(gid);
      const double cn = std::abs(data.Cn()[cn_lid]);

      // compute averaged weighted gap
      const double kappa = cnode->AugData().GetKappa();
      double awgap = cnode->AugData().GetWGap();
      if (kappa != 1.0e12) awgap /= kappa;

      // get normal part of the Lagrange multiplier
      const double zn = cnode->MoData().lm()[0];

      // check nodes of inactive set *************************************
      if (cnode->Active() == false)
      {
        // check for fulfillment of contact condition
        if (zn - cn * awgap > 0.0)
        {
          cnode->Active() = true;
          istatus[ilid] = Status::changed;
        }
      }
      // check nodes of active set ***************************************
      else
      {
        if (zn - cn * awgap <= 0.0)
        {
          cnode->Active() = false;
          istatus[ilid] = Status::changed;
        }
      }
    }
  }  // end loop over all interfaces

  Status lstatus = update_initial_status(cparams, istatus);

  // make local active set status global
  enum Status gstatus = Status::unevaluated;
  {
    int local = static_cast<int>(lstatus);
    int global = static_cast<int>(gstatus);
    strategy_.Comm().MaxAll(&local, &global, 1);
    gstatus = (global > 0 ? Status::changed : Status::unchanged);
  }

  return gstatus;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::ActiveSet::Status CONTACT::Aug::ActiveSet::update_initial_status(
    const CONTACT::ParamsInterface& cparams, const std::vector<enum Status>& istatus) const
{
  static std::vector<std::vector<std::pair<int, bool>>> init_active_list;

  if (not cparams.IsPredictor() or cparams.GetStepNp() != cparams.GetRestartStep() + 1)
  {
    init_active_list.clear();
    return merge(istatus);
  }

  plain_interface_set& interfaces = strategy_.interface_;

  if (init_active_list.size() == 0) init_active_list.resize(interfaces.size());

  std::vector<enum Status> new_istatus(istatus);

  // local interface id
  unsigned ilid = 0;
  for (auto cit = interfaces.cbegin(); cit != interfaces.cend(); ++cit, ++ilid)
  {
    Teuchos::RCP<const CONTACT::Aug::Interface> aug_i_ptr =
        Teuchos::rcp_dynamic_cast<const CONTACT::Aug::Interface>(*cit, true);
    const CONTACT::Aug::Interface& interface = *aug_i_ptr;

    enum class InitStatus
    {
      changed,
      unchanged,
      undefined
    };
    InitStatus initstatus = InitStatus::undefined;
    auto& my_active_node_list = init_active_list[ilid];

    // loop over all slave nodes of the current interface
    const int num_my_slave_row_nodes = interface.SlaveRowNodes()->NumMyElements();
    const int* my_slave_gids = interface.SlaveRowNodes()->MyGlobalElements();

    if (my_active_node_list.size() > 0)
    {
      initstatus = InitStatus::unchanged;
      for (int j = 0; j < num_my_slave_row_nodes; ++j)
      {
        const int sgid = my_slave_gids[j];
        Node* cnode = dynamic_cast<Node*>(interface.Discret().gNode(sgid));
        if (!cnode) FOUR_C_THROW("Cannot find node with gid %", sgid);

        const std::pair<int, bool>& active_pair = my_active_node_list[j];
        if (cnode->Id() != active_pair.first) FOUR_C_THROW("GID mismatch!");
        if (active_pair.second != cnode->Active())
        {
          initstatus = InitStatus::changed;
          break;
        }
      }
    }
    else
    {
      my_active_node_list.resize(num_my_slave_row_nodes);
      for (int j = 0; j < num_my_slave_row_nodes; ++j)
      {
        const int sgid = my_slave_gids[j];
        Node& cnode = static_cast<Node&>(*interface.Discret().gNode(sgid));
        std::pair<int, bool>& active_pair = my_active_node_list[j];
        active_pair = std::make_pair(cnode.Id(), cnode.Active());
      }
    }

    /* Do not set any node initially active, if the default status identification
     * shows a different result compared to the previous run. */
    if (initstatus == InitStatus::changed) continue;

    unsigned set_init_active = 0;
    for (int j = 0; j < num_my_slave_row_nodes; ++j)
    {
      const int sgid = my_slave_gids[j];
      Node& cnode = static_cast<Node&>(*interface.Discret().gNode(sgid));

      // skip nodes which are already active
      if (cnode.Active()) continue;

      /* If this is the first attempt, the node is newly set active but the
       * initial condition and, thus, the interface status must be set to
       * changed. If it is already the 2nd attempt, and the previous test
       * indicates that the active set did not change, the initial condition
       * will lead to exactly the same active set as in the previous execution
       * and, thus, the interface status does not change. */
      if (interface.set_node_initially_active(cparams, cnode) and
          initstatus == InitStatus::undefined)
      {
        new_istatus[ilid] = Status::changed;
        ++set_init_active;
      }
    }

    Core::IO::cout << std::string(60, '*') << "\n"
                   << set_init_active << " slave nodes of interface #" << ilid
                   << " have been set initially active "
                      "via condition line or INITCONTACTBYGAP.\n"
                   << std::string(60, '*') << "\n";
  }

  return merge(new_istatus);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::ActiveSet::Status CONTACT::Aug::ActiveSet::merge(
    const std::vector<Status>& istatus) const
{
  Status status = Status::unevaluated;
  for (const Status is : istatus)
  {
    switch (is)
    {
      case Status::changed:
        return Status::changed;
      case Status::unchanged:
        status = Status::unchanged;
        break;
      case Status::unevaluated:
        FOUR_C_THROW(
            "The active set status is unevaluated! First update the "
            "status before you call this merge routine.");
        exit(EXIT_FAILURE);
      default:
        FOUR_C_THROW("Unknown active set status: %d", is);
        exit(EXIT_FAILURE);
    }
  }

  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::ActiveSet::Print(std::ostream& os) const
{
  plain_interface_set& interfaces = strategy_.interfaces();

  unsigned ilid = 0;
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end(); ++cit)
  {
    os << "---------------- INTERFACE " << ++ilid << " ----------------\n";
    const Interface& interface = dynamic_cast<const Interface&>(**cit);

    // loop over all slave nodes of the current interface
    const int num_my_slave_row_nodes = interface.SlaveRowNodes()->NumMyElements();
    int* my_slave_row_node_gids = interface.SlaveRowNodes()->MyGlobalElements();
    for (int j = 0; j < num_my_slave_row_nodes; ++j)
    {
      const int gid = my_slave_row_node_gids[j];
      Node* cnode = dynamic_cast<Node*>(interface.Discret().gNode(gid));
      if (!cnode) FOUR_C_THROW("Cannot find node with gid %", gid);

      // compute averaged weighted gap
      const double kappa = cnode->AugData().GetKappa();
      double awgap = cnode->AugData().GetWGap();
      if (kappa != 1.0e12) awgap /= kappa;

      // get normal part of the Lagrange multiplier
      const double zn = cnode->MoData().lm()[0];

      os << "PROC #" << strategy_.Comm().MyPID() << " -- NODE #" << std::setprecision(4)
         << std::setw(10) << cnode->Id() << "(" << ilid << ") | " << std::scientific
         << std::setw(12) << cnode->X()[0] << ", " << std::scientific << std::setw(12)
         << cnode->X()[1] << ", " << std::scientific << std::setw(12) << cnode->X()[2]
         << " | wgap = " << std::scientific << std::setw(12) << cnode->AugData().GetWGap()
         << ", awgap = " << std::setw(12) << awgap
         << (cnode->Active() ? "   [ACTIVE]" : " [INACTIVE]") << ", lm = " << std::scientific
         << std::setw(12) << zn << std::endl;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::ActiveSet::update_maps(const CONTACT::ParamsInterface& cparams)
{
  DataContainer& data = strategy_.data();
  plain_interface_set& interfaces = strategy_.interfaces();

  // only if it's a full Newton step...
  if (cparams.IsDefaultStep())
  {
    // store the previous augmented active set
    if (data.g_active_node_row_map_ptr() != Teuchos::null)
      data.g_old_active_slave_nodes_ptr() = Teuchos::rcp(new Epetra_Map(data.GActiveNodeRowMap()));
    else
      data.g_old_active_slave_nodes_ptr() = Teuchos::rcp(new Epetra_Map(0, 0, strategy_.Comm()));
  }
  else
    Core::IO::cout << "This is no default step! History information stays untouched."
                   << Core::IO::endl;

  // (re)setup of the global Epetra_maps
  data.g_active_node_row_map_ptr() = Teuchos::null;
  data.GActiveDofRowMapPtr() = Teuchos::null;
  data.g_active_n_dof_row_map_ptr() = Teuchos::null;
  data.g_active_t_dof_row_map_ptr() = Teuchos::null;

  // loop over all interfaces
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end(); ++cit)
  {
    CONTACT::Interface& interface = **cit;

    // update active set Epetra_Maps
    interface.BuildActiveSet();

    // Update Active set
    data.g_active_node_row_map_ptr() =
        Core::LinAlg::MergeMap(data.g_active_node_row_map_ptr(), interface.ActiveNodes(), false);
    data.GActiveDofRowMapPtr() =
        Core::LinAlg::MergeMap(data.GActiveDofRowMapPtr(), interface.ActiveDofs(), false);
    data.g_active_n_dof_row_map_ptr() =
        Core::LinAlg::MergeMap(data.g_active_n_dof_row_map_ptr(), interface.ActiveNDofs(), false);
    data.g_active_t_dof_row_map_ptr() =
        Core::LinAlg::MergeMap(data.g_active_t_dof_row_map_ptr(), interface.ActiveTDofs(), false);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::ActiveSet::sanity_check(
    const CONTACT::ParamsInterface& cparams, const enum Status gstatus) const
{
  const DataContainer& data = strategy_.data();

  if (cparams.IsDefaultStep() and (gstatus == Status::changed) == data.is_active_set_converged())
    FOUR_C_THROW(
        "The convergence state of the active set has not been correctly "
        "detected: %s",
        status2_string(gstatus).c_str());

  Core::IO::cout << "old number of active nodes:     "
                 << data.g_old_active_slave_nodes_ptr()->NumGlobalElements() << Core::IO::endl;
  Core::IO::cout << "current number of active nodes: "
                 << data.g_active_node_row_map_ptr()->NumGlobalElements() << Core::IO::endl;
}

FOUR_C_NAMESPACE_CLOSE
