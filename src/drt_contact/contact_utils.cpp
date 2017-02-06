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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::UTILS::GetContactConditions(
    std::vector<DRT::Condition*>& ccond,
    const std::vector<DRT::Condition*> beamandsolidcontactconditions,
    const bool & throw_error)
{
  /* Sort out beam-to-solid contact pairs, since these are treated in the
   * beam3contact framework */
  for (std::size_t i = 0; i < beamandsolidcontactconditions.size(); ++i)
  {
    if(*(beamandsolidcontactconditions[i]->Get<std::string>("Application"))
        !="Beamtosolidcontact")
    {
      ccond.push_back(beamandsolidcontactconditions[i]);
    }
  }

  /* There must be more than one contact condition unless we have a self
   * contact problem! */
  if (ccond.size() < 1)
  {
    if (throw_error)
      dserror("ERROR: Not enough contact conditions in discretization");
    return -1;
  }
  if (ccond.size() == 1)
  {
    const std::string* side = ccond[0]->Get<std::string>("Side");
    if (*side != "Selfcontact")
    {
      if (throw_error)
        dserror("ERROR: Not enough contact conditions in discretization");
      return -2;
    }
  }
  // everything worked fine
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::UTILS::GetContactConditionGroups(
    std::vector<std::vector<DRT::Condition*> >& ccond_grps,
    const DRT::DiscretizationInterface& discret,
    const bool & throw_error)
{
  // vector that contains solid-to-solid and beam-to-solid contact pairs
  std::vector<DRT::Condition*> beamandsolidcontactconditions(0);
  discret.GetCondition("Contact", beamandsolidcontactconditions);

  std::vector<DRT::Condition*> cconds(0);
  int err = CONTACT::UTILS::GetContactConditions(cconds,
      beamandsolidcontactconditions,throw_error);
  // direct return, if an error occurred
  if (err)
    return err;
  CONTACT::UTILS::GetContactConditionGroups(ccond_grps,cconds);
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::UTILS::GetContactConditionGroups(
    std::vector<std::vector<DRT::Condition*> >& ccond_grps,
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
    const std::vector<int>* group1v = current_grp[0]->Get<std::vector<int> >(
        "Interface ID");
    if (!group1v)
      dserror("ERROR: Contact Conditions does not have value 'Interface ID'");
    int groupid1 = (*group1v)[0];
    bool foundit = false;

    // only one surface per group is ok for self contact
    const std::string* side = cconds[i]->Get<std::string>("Side");
    if (*side == "Selfcontact")
      foundit = true;

    for (std::size_t j = 0; j < cconds.size(); ++j)
    {
      // do not compare ids of one and the same contact condition
      if (j == i)
        continue;
      tempcond = cconds[j];
      const std::vector<int>* group2v = tempcond->Get<std::vector<int> >(
          "Interface ID");
      if (!group2v)
        dserror("ERROR: Contact Conditions does not have value 'Interface ID'");
      int groupid2 = (*group2v)[0];
      // Do the IDs coincide?
      if (groupid1 != groupid2)
        continue; // not in the group
      foundit = true; // found a group entry
      current_grp.push_back(tempcond); // store it in the current group
    }

    // now we should have found a group of conditions
    if (!foundit)
      dserror("ERROR: Cannot find matching contact condition for id %i", groupid1);

    // see whether we found this group before
    bool foundbefore = false;
    for (std::size_t j = 0; j < found_grps.size(); ++j)
      if (groupid1 == found_grps[j])
      {
        foundbefore = true;
        break;
      }

    // if we have processed this group before, do nothing
    if (foundbefore)
      continue;

    // we have not found this group before, store it
    found_grps.push_back(groupid1);

    // store the new unique group of contact conditions
    ccond_grps.push_back(current_grp);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::UTILS::GetMasterSlaveSideInfo(
    std::vector<bool>& isslave,
    std::vector<bool>& isself,
    const std::vector<DRT::Condition*> cond_grp)
{
  bool hasslave = false;
  bool hasmaster = false;
  std::vector<const std::string*> sides(cond_grp.size());
  // safety...
  isslave.clear();
  isslave.resize(cond_grp.size(),false);
  isself.clear();
  isself.resize(cond_grp.size(),false);

  for (int j = 0; j < (int) sides.size(); ++j)
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

  if (!hasslave)
    dserror("ERROR: Slave side missing in contact condition group!");
  if (!hasmaster)
    dserror("ERROR: Master side missing in contact condition group!");

  // check for self contact group
  if (isself[0])
  {
    for (int j = 1; j < (int) isself.size(); ++j)
      if (!isself[j])
        dserror("ERROR: Inconsistent definition of self contact condition group!");
  }
}
