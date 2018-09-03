/*----------------------------------------------------------------------------*/
/**
\file contact_utils.cpp

\brief Contains a summary of contact utility functions

\maintainer Michael Hiermeier

\date Jun 15, 2016

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "contact_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_io/every_iteration_writer.H"

#include "../drt_inpar/inpar_mortar.H"

#include "../linalg/linalg_serialdensematrix.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::UTILS::GetContactConditions(std::vector<DRT::Condition*>& ccond,
    const std::vector<DRT::Condition*> beamandsolidcontactconditions, const bool& throw_error)
{
  /* Sort out beam-to-solid contact pairs, since these are treated in the
   * beam3contact framework */
  for (std::size_t i = 0; i < beamandsolidcontactconditions.size(); ++i)
  {
    if (*(beamandsolidcontactconditions[i]->Get<std::string>("Application")) !=
        "Beamtosolidcontact")
    {
      ccond.push_back(beamandsolidcontactconditions[i]);
    }
  }

  /* There must be more than one contact condition unless we have a self
   * contact problem! */
  if (ccond.size() < 1)
  {
    if (throw_error) dserror("ERROR: Not enough contact conditions in discretization");
    return -1;
  }
  if (ccond.size() == 1)
  {
    const std::string* side = ccond[0]->Get<std::string>("Side");
    if (*side != "Selfcontact")
    {
      if (throw_error) dserror("ERROR: Not enough contact conditions in discretization");
      return -2;
    }
  }
  // everything worked fine
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::UTILS::GetContactConditionGroups(std::vector<std::vector<DRT::Condition*>>& ccond_grps,
    const DRT::DiscretizationInterface& discret, const bool& throw_error)
{
  // vector that contains solid-to-solid and beam-to-solid contact pairs
  std::vector<DRT::Condition*> beamandsolidcontactconditions(0);
  discret.GetCondition("Contact", beamandsolidcontactconditions);

  std::vector<DRT::Condition*> cconds(0);
  int err =
      CONTACT::UTILS::GetContactConditions(cconds, beamandsolidcontactconditions, throw_error);
  // direct return, if an error occurred
  if (err) return err;
  CONTACT::UTILS::GetContactConditionGroups(ccond_grps, cconds);
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::UTILS::GetContactConditionGroups(
    std::vector<std::vector<DRT::Condition*>>& ccond_grps,
    const std::vector<DRT::Condition*>& cconds)
{
  ccond_grps.clear();
  /* find all pairs of matching contact conditions
   * there is a maximum of (conditions / 2) groups */
  std::vector<int> found_grps(0);

  for (std::size_t i = 0; i < cconds.size(); ++i)
  {
    std::vector<DRT::Condition*> current_grp(0);
    DRT::Condition* tempcond = NULL;

    // try to build contact group around this condition
    current_grp.push_back(cconds[i]);
    const std::vector<int>* group1v = current_grp[0]->Get<std::vector<int>>("Interface ID");
    if (!group1v) dserror("ERROR: Contact Conditions does not have value 'Interface ID'");
    int groupid1 = (*group1v)[0];
    bool foundit = false;

    // only one surface per group is ok for self contact
    const std::string* side = cconds[i]->Get<std::string>("Side");
    if (*side == "Selfcontact") foundit = true;

    for (std::size_t j = 0; j < cconds.size(); ++j)
    {
      // do not compare ids of one and the same contact condition
      if (j == i) continue;
      tempcond = cconds[j];
      const std::vector<int>* group2v = tempcond->Get<std::vector<int>>("Interface ID");
      if (!group2v) dserror("ERROR: Contact Conditions does not have value 'Interface ID'");
      int groupid2 = (*group2v)[0];
      // Do the IDs coincide?
      if (groupid1 != groupid2) continue;  // not in the group
      foundit = true;                      // found a group entry
      current_grp.push_back(tempcond);     // store it in the current group
    }

    // now we should have found a group of conditions
    if (!foundit) dserror("ERROR: Cannot find matching contact condition for id %i", groupid1);

    // see whether we found this group before
    bool foundbefore = false;
    for (std::size_t j = 0; j < found_grps.size(); ++j)
      if (groupid1 == found_grps[j])
      {
        foundbefore = true;
        break;
      }

    // if we have processed this group before, do nothing
    if (foundbefore) continue;

    // we have not found this group before, store it
    found_grps.push_back(groupid1);

    // store the new unique group of contact conditions
    ccond_grps.push_back(current_grp);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::UTILS::GetMasterSlaveSideInfo(std::vector<bool>& isslave, std::vector<bool>& isself,
    const std::vector<DRT::Condition*> cond_grp)
{
  bool hasslave = false;
  bool hasmaster = false;
  std::vector<const std::string*> sides(cond_grp.size());
  // safety...
  isslave.clear();
  isslave.resize(cond_grp.size(), false);
  isself.clear();
  isself.resize(cond_grp.size(), false);

  for (int j = 0; j < (int)sides.size(); ++j)
  {
    sides[j] = cond_grp[j]->Get<std::string>("Side");
    if (*sides[j] == "Slave")
    {
      hasslave = true;
      isslave[j] = true;
      isself[j] = false;
    }
    else if (*sides[j] == "Master")
    {
      hasmaster = true;
      isslave[j] = false;
      isself[j] = false;
    }
    else if (*sides[j] == "Selfcontact")
    {
      hasmaster = true;
      hasslave = true;
      isslave[j] = false;
      isself[j] = true;
    }
    else
      dserror("ERROR: Unknown contact side qualifier!");
  }

  if (!hasslave) dserror("ERROR: Slave side missing in contact condition group!");
  if (!hasmaster) dserror("ERROR: Master side missing in contact condition group!");

  // check for self contact group
  if (isself[0])
  {
    for (int j = 1; j < (int)isself.size(); ++j)
      if (!isself[j]) dserror("ERROR: Inconsistent definition of self contact condition group!");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::UTILS::WriteConservationDataToFile(const int mypid, const int interface_id,
    const int nln_iter, const LINALG::SerialDenseMatrix& conservation_data,
    const std::string& ofile_path, const std::string& prefix)
{
  if (mypid != 0) return;

  static std::vector<std::string> done_prefixes;

  const std::string path(IO::ExtractPath(ofile_path));
  const std::string dir_name(
      IO::RemoveRestartCountFromFileName(IO::ExtractFileName(ofile_path)) + "_conservation");

  std::string full_filepath(path + dir_name);
  IO::CreateDirectory(full_filepath, mypid);
  full_filepath += "/" + prefix + "_" + "conservation.data";

  bool is_done = false;
  for (const std::string& done_prefix : done_prefixes)
  {
    if (done_prefix == prefix)
    {
      is_done = true;
      break;
    }
  }

  // first attempt: clear file content and write header
  if (not is_done)
  {
    done_prefixes.push_back(prefix);

    std::ofstream of(full_filepath, std::ios_base::out);
    of << std::setw(24) << "it" << std::setw(24) << "interface" << std::setw(24) << "Fsl_X"
       << std::setw(24) << "Fsl_Y" << std::setw(24) << "Fsl_Z" << std::setw(24) << "Fma_X"
       << std::setw(24) << "Fma_Y" << std::setw(24) << "Fma_Z" << std::setw(24) << "Fb_X"
       << std::setw(24) << "Fb_Y" << std::setw(24) << "Fb_Z" << std::setw(24) << "Mosl_X"
       << std::setw(24) << "Mosl_Y" << std::setw(24) << "Mosl_Z" << std::setw(24) << "Moma_X"
       << std::setw(24) << "Moma_Y" << std::setw(24) << "Moma_Z" << std::setw(24) << "Mob_X"
       << std::setw(24) << "Mob_Y" << std::setw(24) << "Mob_Z\n";
    of.close();
  }

  std::ofstream of(full_filepath, std::ios_base::out | std::ios_base::app);

  if (conservation_data.M() < 18) dserror("The conservation_data has insufficient size!");

  of << std::setw(24) << nln_iter << std::setw(24) << interface_id;
  of << std::setprecision(16);
  for (unsigned i = 0; i < static_cast<unsigned>(conservation_data.M()); ++i)
  {
    of << std::setw(24) << std::setw(24) << std::scientific << conservation_data(i, 0);
  }
  of << "\n";
  of.close();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::UTILS::DbcHandler::DetectDbcSlaveNodesAndElements(
    const DRT::DiscretizationInterface& str_discret,
    const std::vector<std::vector<DRT::Condition*>>& ccond_grps,
    std::set<const DRT::Node*>& dbc_slave_nodes, std::set<const DRT::Element*>& dbc_slave_eles)
{
  dbc_slave_nodes.clear();
  dbc_slave_eles.clear();

  std::map<const DRT::Node*, int> dbc_slave_node_map;

  std::vector<const DRT::Condition*> sl_conds;

  for (auto& ccond_grp : ccond_grps)
  {
    std::vector<bool> isslave;
    std::vector<bool> isself;
    CONTACT::UTILS::GetMasterSlaveSideInfo(isslave, isself, ccond_grp);

    for (unsigned i = 0; i < ccond_grp.size(); ++i)
    {
      if (not isslave[i]) continue;

      const DRT::Condition* sl_cond = ccond_grp[i];

      const int dbc_handling_id = sl_cond->GetInt("dbc_handling");
      const INPAR::MORTAR::DBCHandling dbc_handling =
          static_cast<INPAR::MORTAR::DBCHandling>(dbc_handling_id);

      switch (dbc_handling)
      {
        case INPAR::MORTAR::DBCHandling::remove_dbc_nodes_from_slave_side:
        {
          sl_conds.push_back(sl_cond);
          break;
        }
        case INPAR::MORTAR::DBCHandling::do_nothing:
        {
          break;
        }
        default:
          dserror("Unknown dbc_handlin enum %d", dbc_handling);
          exit(EXIT_FAILURE);
      }
    }
  }

  DetectDbcSlaveNodes(dbc_slave_node_map, str_discret, sl_conds);

  for (auto& dbc_slave_node_pair : dbc_slave_node_map)
    dbc_slave_nodes.insert(dbc_slave_node_pair.first);

  DetectDbcSlaveElements(dbc_slave_eles, dbc_slave_node_map, sl_conds);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::UTILS::DbcHandler::DetectDbcSlaveNodes(
    std::map<const DRT::Node*, int>& dbc_slave_node_map,
    const DRT::DiscretizationInterface& str_discret,
    const std::vector<const DRT::Condition*>& sl_conds)
{
  std::vector<DRT::Condition*> dconds;
  str_discret.GetCondition("Dirichlet", dconds);

  // collect all slave node ids
  std::vector<std::pair<int, int>> slnodeids;
  for (auto& sl_cond : sl_conds)
  {
    const std::vector<int>* sl_nids = sl_cond->Get<std::vector<int>>("Node Ids");
    slnodeids.reserve(slnodeids.size() + sl_nids->size());
    for (int sl_nid : *sl_nids) slnodeids.push_back(std::make_pair(sl_nid, sl_cond->Id()));
  }

  for (std::pair<int, int> slpair : slnodeids)
  {
    const int snid = slpair.first;

    bool found = false;

    for (DRT::Condition* dcond : dconds)
    {
      const std::vector<int>* dnids = dcond->Get<std::vector<int>>("Node Ids");
      for (int dnid : *dnids)
      {
        if (snid == dnid)
        {
          found = true;
          break;
        }
      }
      if (found) break;
    }

    // skip non dbc nodes
    if (not found) continue;

    // skip NULL ptrs
    if (str_discret.HaveGlobalNode(snid))
    {
      const DRT::Node* node = str_discret.gNode(snid);
      dbc_slave_node_map.insert(std::make_pair(node, slpair.second));
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::UTILS::DbcHandler::DetectDbcSlaveElements(
    std::set<const DRT::Element*>& dbc_slave_eles,
    const std::map<const DRT::Node*, int>& dbc_slave_nodes,
    const std::vector<const DRT::Condition*>& sl_conds)
{
  for (auto& dbc_sl_node : dbc_slave_nodes)
  {
    const int slnid = dbc_sl_node.first->Id();
    const int slcond_id = dbc_sl_node.second;

    auto sl_citer = sl_conds.cbegin();
    while (sl_citer != sl_conds.cend())
    {
      if ((*sl_citer)->Id() == slcond_id) break;

      ++sl_citer;
    }
    const DRT::Condition& slcond = **sl_citer;

    const std::map<int, Teuchos::RCP<DRT::Element>>& geometry = slcond.Geometry();
    for (auto& iele_pair : geometry)
    {
      const DRT::Element* ele = iele_pair.second.get();

      const int* ele_nids = ele->NodeIds();
      for (int i = 0; i < ele->NumNode(); ++i)
      {
        const int ele_nid = ele_nids[i];

        if (ele_nid == slnid) dbc_slave_eles.insert(ele);
      }
    }
  }
}
