/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "ssi_utils.H"

#include <Epetra_Map.h>

#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_adapter/adapter_coupling.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 09/2014 */
/* Function for checking that the different time steps are a
 multiplicative of each other                                           */

int SSI::Utils::CheckTimeStepping(double dt1, double dt2)
{
  double workdt1 = std::min(dt1, dt2);
  double workdt2 = std::max(dt1, dt2);
  double t1 = 0.0;
  int i = 0;

  while (true)
  {
    i++;
    t1 = i * workdt1;

    if (std::abs(t1 - workdt2) < 10E-10)
      break;

    else if (t1 > workdt2)
      dserror("Chosen time steps %f and %f are not a multiplicative of each other", dt1, dt2);
  }
  return i;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 10/2014 */
// Modification of time parameter list for problem with different time step size

void SSI::Utils::ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ssiparams,
    Teuchos::ParameterList& scatradyn, Teuchos::ParameterList& sdyn)
{
  bool difftimestep = DRT::INPUT::IntegralValue<int>(ssiparams, "DIFFTIMESTEPSIZE");

  if (difftimestep)  // Create subproblems with different time steps
  {
    // Check correct choice of time stepping for single fields
    double scatrastep = scatradyn.get<double>("TIMESTEP");
    double solidstep = sdyn.get<double>("TIMESTEP");

    SSI::Utils::CheckTimeStepping(scatrastep, solidstep);

    // modify global time step size
    ssiparams.set<double>("TIMESTEP", std::min(scatrastep, solidstep));
  }
  else
  {
    // -------------------------------------------------------------------
    // overrule certain parameters for coupled problems
    // -------------------------------------------------------------------
    // the default time step size
    scatradyn.set<double>("TIMESTEP", ssiparams.get<double>("TIMESTEP"));
    sdyn.set<double>("TIMESTEP", ssiparams.get<double>("TIMESTEP"));
    // maximum simulation time
    scatradyn.set<double>("MAXTIME", ssiparams.get<double>("MAXTIME"));
    sdyn.set<double>("MAXTIME", ssiparams.get<double>("MAXTIME"));
    // maximum number of timesteps
    scatradyn.set<int>("NUMSTEP", ssiparams.get<int>("NUMSTEP"));
    sdyn.set<int>("NUMSTEP", ssiparams.get<int>("NUMSTEP"));
  }

  // Check correct input of restart. Code relies that both time value RESTARTEVRYTIME and
  // RESULTSEVRYTIME are given if restart from time is applied
  double restarttime = ssiparams.get<double>("RESTARTEVRYTIME");
  double updatetime = ssiparams.get<double>("RESULTSEVRYTIME");
  if ((updatetime > 0.0) or (restarttime > 0.0))
    if (!(updatetime > 0.0) and !(restarttime > 0.0))
      dserror(
          "If time controlled output and restart is desired, both parameters RESTARTEVRYTIME and "
          "RESULTSEVRYTIME has to be set");

  // set restart params
  int scatrarestart;
  int structurerestart;

  if (restarttime > 0.0)
  {
    scatrarestart = SSI::Utils::CheckTimeStepping(scatradyn.get<double>("TIMESTEP"), restarttime);
    structurerestart = SSI::Utils::CheckTimeStepping(sdyn.get<double>("TIMESTEP"), restarttime);
  }
  else
  {
    int restart = ssiparams.get<int>("RESTARTEVRY");
    scatrarestart = restart;
    structurerestart = restart;
  }

  // set output params
  int scatraupres;
  int structureupres;

  if (updatetime > 0.0)
  {
    scatraupres = SSI::Utils::CheckTimeStepping(scatradyn.get<double>("TIMESTEP"), updatetime);
    structureupres = SSI::Utils::CheckTimeStepping(sdyn.get<double>("TIMESTEP"), updatetime);
  }
  else
  {
    int update = ssiparams.get<int>("RESULTSEVRY");
    scatraupres = update;
    structureupres = update;
  }

  // restart
  scatradyn.set<int>("RESTARTEVRY", scatrarestart);
  sdyn.set<int>("RESTARTEVRY", structurerestart);
  // solution output
  scatradyn.set<int>("RESULTSEVRY", scatraupres);
  sdyn.set<int>("RESULTSEVRY", structureupres);

  if (comm.MyPID() == 0)
  {
    std::cout << "====================== Overview of chosen time stepping: "
                 "==============================\n"
              << "\t Timestep scatra:           " << scatradyn.get<double>("TIMESTEP") << "\n"
              << "\t Timestep structure:        " << sdyn.get<double>("TIMESTEP") << "\n"
              << "\t Result step scatra:        " << scatradyn.get<int>("RESULTSEVRY") << "\n"
              << "\t Result step structure:     " << sdyn.get<int>("RESULTSEVRY") << "\n"
              << "\t Restart step scatra:       " << scatradyn.get<int>("RESTARTEVRY") << "\n"
              << "\t Restart step structure:    " << sdyn.get<int>("RESTARTEVRY") << "\n"
              << "================================================================================="
                 "=======\n \n";
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::Utils::CheckConsistencyWithS2IMeshtyingCondition(
    std::vector<DRT::Condition*> conditionsToBeTested, Teuchos::RCP<DRT::Discretization>& structdis)
{
  std::vector<DRT::Condition*> s2iconditions;
  structdis->GetCondition("S2ICoupling", s2iconditions);

  // loop over all conditions to be tested and check for a consistent initialization of the s2i
  // conditions
  for (auto conditionToBeTested : conditionsToBeTested)
  {
    bool matchingconditions(false);
    bool isslave(true);
    const int s2icouplingid = conditionToBeTested->GetInt("S2ICouplingID");
    auto* side = conditionToBeTested->Get<std::string>("Side");
    // check interface side
    if (*side == "Slave")
      isslave = true;
    else if (*side == "Master")
      isslave = false;
    else
      dserror(
          "Interface side of tested condition not recognized, has to be either 'Slave' or "
          "'Master'");

    // loop over all s2i conditions to find the one that is matching the current ssi condition
    for (auto s2icondition : s2iconditions)
    {
      const int s2iconditionid = s2icondition->GetInt("ConditionID");
      // only do further checks if Ids match
      if (s2icouplingid != s2iconditionid) continue;

      // check the interface side
      switch (s2icondition->GetInt("interface side"))
      {
        case INPAR::S2I::side_slave:
        {
          if (isslave)
            matchingconditions = DRT::UTILS::HaveSameNodes(conditionToBeTested, s2icondition);

          break;
        }
        case INPAR::S2I::side_master:
        {
          if (!isslave)
            matchingconditions = DRT::UTILS::HaveSameNodes(conditionToBeTested, s2icondition);

          break;
        }
        default:
        {
          dserror("interface side of 'S2iCondition' has to be either 'Slave' or 'Master'");
          break;
        }
      }
    }

    if (!matchingconditions)
      dserror(
          "Did not find 'S2ICoupling' condition with ID: %i and interface side: %s as defined in "
          "the condition to be tested",
          s2icouplingid, side->c_str());
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::Coupling> SSI::Utils::SetupInterfaceCouplingAdapterStructure(
    Teuchos::RCP<DRT::Discretization> structdis)
{
  // initialize integer vectors for global IDs of master-side and slave-side interface nodes on
  // structure discretization
  std::vector<int> inodegidvec_master;
  std::vector<int> inodegidvec_slave;

  // extract scatra-scatra interface coupling conditions from structure discretization
  std::vector<DRT::Condition*> conditions(0, nullptr);
  structdis->GetCondition("S2ICoupling", conditions);

  // loop over all conditions
  for (const auto& condition : conditions)
  {
    // extract interface side associated with current condition
    const int side = condition->GetInt("interface side");

    // extract nodes associated with current condition
    const std::vector<int>* const inodegids = condition->Nodes();

    for (const int inodegid : *inodegids)
    {
      // insert global ID of current node into associated vector only if node is owned by current
      // processor need to make sure that node is stored on current processor, otherwise cannot
      // resolve "->Owner()"
      if (structdis->HaveGlobalNode(inodegid) and
          structdis->gNode(inodegid)->Owner() == structdis->Comm().MyPID())
        side == INPAR::S2I::side_master ? inodegidvec_master.push_back(inodegid)
                                        : inodegidvec_slave.push_back(inodegid);
    }
  }

  // remove potential duplicates from vectors
  std::sort(inodegidvec_master.begin(), inodegidvec_master.end());
  inodegidvec_master.erase(
      unique(inodegidvec_master.begin(), inodegidvec_master.end()), inodegidvec_master.end());
  std::sort(inodegidvec_slave.begin(), inodegidvec_slave.end());
  inodegidvec_slave.erase(
      unique(inodegidvec_slave.begin(), inodegidvec_slave.end()), inodegidvec_slave.end());

  // setup scatra-scatra interface coupling adapter for structure field
  auto coupling_structure = Teuchos::rcp(new ADAPTER::Coupling());
  coupling_structure->SetupCoupling(*structdis, *structdis, inodegidvec_master, inodegidvec_slave,
      DRT::Problem::Instance()->NDim(), true, 1.0e-8);

  return coupling_structure;
}
