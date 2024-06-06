/*----------------------------------------------------------------------*/
/*! \file
  \brief  Utility Methods For Fluid Porous Structure Interaction Problems
\level 3

 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  rauch 11/12 |
 *----------------------------------------------------------------------*/

// FPSI includes
#include "4C_fpsi_utils.hpp"

#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_discretization_condition_selector.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_fpsi_monolithic_plain.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fpsi.hpp"
#include "4C_poroelast_scatra_utils_clonestrategy.hpp"
#include "4C_poroelast_scatra_utils_setup.hpp"
#include "4C_poroelast_utils_clonestrategy.hpp"
#include "4C_poroelast_utils_setup.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

Teuchos::RCP<FPSI::Utils> FPSI::Utils::instance_;

/*----------------------------------------------------------------------/
| Instance method for singleton class                            rauch  |
/----------------------------------------------------------------------*/
Teuchos::RCP<FPSI::Utils> FPSI::Utils::Instance()
{
  if (instance_ == Teuchos::null) instance_ = Teuchos::rcp(new Utils());

  return instance_;
}

/*----------------------------------------------------------------------/
| Setup discretization                                           rauch  |
/----------------------------------------------------------------------*/
Teuchos::RCP<FPSI::FpsiBase> FPSI::Utils::setup_discretizations(const Epetra_Comm& comm,
    const Teuchos::ParameterList& fpsidynparams, const Teuchos::ParameterList& poroelastdynparams)
{
  Global::Problem* problem = Global::Problem::Instance();

  fluid_poro_fluid_interface_map_ = Teuchos::rcp(new std::map<int, int>);
  poro_fluid_fluid_interface_map_ = Teuchos::rcp(new std::map<int, int>);

  // 1.-Initialization.
  Teuchos::RCP<Discret::Discretization> structdis = problem->GetDis("structure");
  if (structdis == Teuchos::null) FOUR_C_THROW(" !!! structdis empty !!! Awwww MAAAAN !!!");
  Teuchos::RCP<Discret::Discretization> porofluiddis = problem->GetDis("porofluid");
  if (porofluiddis == Teuchos::null) FOUR_C_THROW(" !!! porofluiddis empty !!! Awwww MAAAAN !!!");
  Teuchos::RCP<Discret::Discretization> fluiddis = problem->GetDis("fluid");
  if (fluiddis == Teuchos::null) FOUR_C_THROW(" !!! fluiddis empty !!! Awwww MAAAAN !!!");
  Teuchos::RCP<Discret::Discretization> aledis = problem->GetDis("ale");
  if (aledis == Teuchos::null) FOUR_C_THROW(" !!! aledis empty !!! Awwww MAAAAN !!!");

  /*
   Ok Dude, listen carefully for the following is very important. The stuff below is
   actually quite basic but the emphasis is put on the right order of doing it. In
   general the filling order of the discretizations determines their
   individual dof numbers and thus the mapping of the local indices to the global ones.
   The rule is simple: First filled first mapped. So have a look at what we do below.
   First we fill the structural discretization which means that its local indices equal
   its global indices beginning from 0 and ending with the number of structural dofs
   depending on the specific problem size under consideration. Now the porofluiddis is
   filled. Note that this discretization doesn't exist at that moment because in our
   input file we only have a structural discretization for the domain of the porous
   medium ("structure") and a fluid discretization for the domain containing the free
   fluid ("fluid"). To bring the porofluiddis to life we clone it from the structdis.
   Since it was filled right after the structdis the dof numbering of the porous structure
   and the fluid contained by it are adjacent. Now we fill the fluiddis (which exists) and
   the aledis (which is nonexistent just like the porofluiddis before the cloning procedure).
   We expect the fluiddis to have smaller dof numbers than the ale discretization. Finally
   cloning the aledis from the fluiddis makes this dream come true. Building the monolithic
   systemmatrix based on the DofMaps built before, results in the following order of the
   diagonal blocks (from northwest to southeast): porostructure - porofluid - free fluid - ale.

   Thats nice....

   This would be corrupted however if we were doing the following:

   structdis    -> fill_complete()
   porofluiddis -> fill_complete()
   fluiddis     -> fill_complete()
   aledis       -> fill_complete()

   Clone(structdis,porofluiddis)
   Clone(fluiddis,aledis)

   Thats not so nice because then the fluiddis would have the same DofMap like the porofluiddis.
   The reason for this awkward mess is that when filling the fluiddis after the porofluiddis
   the DofNumbering of the structuredis and the fluiddis will be adjacent for there exists no
   porofluiddis to append the fluiddofmap to. The resulting monolithic blockmatrix will assemble
   into the same positions like the porofluiddis. Consequently not even the matrix dimension
   will be correct. That's because when the porofluiddis is cloned from the structdis in the end the
   DofMap of the porofluiddis will also be adjacent to the structdis DofMap.

   You get it??

   Thats why it's important to immediately clone the nonexistent discretizations after filling
   them.

   */
  // setup of the discretizations, including clone strategy

  // choose cloning strategy depending on poroelast or scatra poroelast problem type
  if (problem->GetProblemType() == Core::ProblemType::fps3i)
  {
    PoroElast::UTILS::SetupPoro<PoroElastScaTra::UTILS::PoroelastCloneStrategyforScatraElements>();
  }
  else
  {
    PoroElast::UTILS::SetupPoro<PoroElast::UTILS::PoroelastCloneStrategy>();
  }


  fluiddis->fill_complete(true, true, true);
  aledis->fill_complete(true, true, true);

  // 3.- Create ALE elements if the ale discretization is empty
  if (aledis->NumGlobalNodes() == 0)  // ALE discretization still empty
  {
    Core::FE::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(
        fluiddis, aledis, Global::Problem::Instance()->CloningMaterialMap());
    aledis->fill_complete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else  // ALE discretization already filled
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis, aledis))
      FOUR_C_THROW(
          "Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use the ALE cloning functionality or ensure non-overlapping node numbering!");
  }

  SetupInterfaceMap(comm, structdis, porofluiddis, fluiddis, aledis);

  // 4.- get coupling algorithm
  Teuchos::RCP<FPSI::FpsiBase> fpsi_algo = Teuchos::null;
  int coupling = Core::UTILS::IntegralValue<int>(fpsidynparams, "COUPALGO");
  switch (coupling)
  {
    case fpsi_monolithic_plain:
    {
      fpsi_algo = Teuchos::rcp(new FPSI::MonolithicPlain(comm, fpsidynparams, poroelastdynparams));
      break;
    }  // case monolithic
    case partitioned:
    {
      FOUR_C_THROW(
          "Partitioned solution scheme not implemented for FPSI, yet. "
          "Make sure that the parameter COUPALGO is set to 'fpsi_monolithic_plain', "
          "and the parameter PARITIONED is set to 'monolithic'. ");
      Inpar::FPSI::PartitionedCouplingMethod method;
      method = Core::UTILS::IntegralValue<Inpar::FPSI::PartitionedCouplingMethod>(
          fpsidynparams, "PARTITIONED");
      if (method == Inpar::FPSI::RobinNeumann)
      {
        // do nothing
      }
      else
      {
        FOUR_C_THROW(
            "Only RobinNeumann algorithm available for partitioned FPSI !!!\n"
            "Set 'PARTITIONED' to 'RobinNeumann' in input file.");
      }
      break;
    }
  }  // switch(coupling)

  return fpsi_algo;
}  // SetupDiscretization()


/*---------------------------------------------------------------------------/
| Setup Local Interface Facing Element Map (for parallel distr.)      rauch  |
/---------------------------------------------------------------------------*/
void FPSI::Utils::setup_local_interface_facing_element_map(Discret::Discretization& masterdis,
    const Discret::Discretization& slavedis, const std::string& condname,
    std::map<int, int>& interfacefacingelementmap)
{
  Global::Problem* problem = Global::Problem::Instance();
  const Epetra_Comm& mastercomm = problem->GetDis(masterdis.Name())->Comm();

  bool condition_exists = true;

  Core::Conditions::Condition* slavecond = slavedis.GetCondition(condname);
  if (slavecond == nullptr)
  {
    condition_exists = false;
    std::cout << "WARNING: Condition <" << condname << "> does not exist on discretisation <"
              << slavedis.Name() << ">! (OK if no " << condname << "-Interface is wanted)"
              << std::endl;
  }

  Core::Conditions::Condition* mastercond = masterdis.GetCondition(condname);
  if (mastercond == nullptr)
  {
    condition_exists = false;
    std::cout << "WARNING: Condition <" << condname << "> does not exist on discretisation <"
              << masterdis.Name() << ">! (OK if no " << condname << "-Interface is wanted)"
              << std::endl;
  }

  if (!condition_exists) return;

  std::map<int, Teuchos::RCP<Core::Elements::Element>>& slavegeom = slavecond->Geometry();
  std::map<int, Teuchos::RCP<Core::Elements::Element>>& mastergeom = mastercond->Geometry();

  std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator curr;
  std::multimap<int, double> slaveinterfaceelementidentificationmap;

  ///////////////////////////////////////////////////////////////////////////////////
  ////////////                                                           ////////////
  ////////////                         SLAVE                             ////////////
  ////////////                       Interface                           ////////////
  ///////////////////////////////////////////////////////////////////////////////////
  int slavegeomsize = slavegeom.size();
  int globalslavegeomsize;
  mastercomm.SumAll(&slavegeomsize, &globalslavegeomsize, 1);

  // do for every slave interface element
  for (curr = slavegeom.begin(); curr != slavegeom.end(); ++curr)
  {
    // std::cout<<"Proc "<<mastercomm.MyPID()<<" owns slavegeom "<<curr->first<<endl;

    std::vector<double> slaveloc;
    slaveloc.assign(3, 0.0);

    //       std::cout<<"Current Slave Interface Element ID: "<<curr->second->Id()<<endl;

    if (slavedis.HaveDofs() == false)
    {
      FOUR_C_THROW("No degrees of freedom have been assigned to discretization");
    }

    // do for every interface slave node
    for (int nodenum = 0; nodenum < curr->second->num_node(); nodenum++)
    {
      const Core::Nodes::Node* const* currslavenode = curr->second->Nodes();

      std::vector<double> temploc;
      temploc.assign(3, 0.0);
      temploc[0] = currslavenode[nodenum]->X()[0];
      temploc[1] = currslavenode[nodenum]->X()[1];
      temploc[2] = currslavenode[nodenum]->X()[2];
      for (int dim = 0; dim < 3; dim++)
      {
        slaveloc[dim] = slaveloc[dim] + temploc[dim];
      }
    }  // for every slave node

    slaveinterfaceelementidentificationmap.insert(
        std::pair<int, double>(curr->second->Id(), slaveloc[0]));
    slaveinterfaceelementidentificationmap.insert(
        std::pair<int, double>(curr->second->Id(), slaveloc[1]));
    slaveinterfaceelementidentificationmap.insert(
        std::pair<int, double>(curr->second->Id(), slaveloc[2]));

  }  // do for every global slave element

  mastercomm.Barrier();  // wait for procs

  //  // printout
  //  for(int proc = 0; proc < (mastercomm.NumProc()+1); proc++)
  //  {
  //    if(mastercomm.MyPID() == proc)
  //    {
  //      std::cout<<"\n slaveinterfaceelementidentificationmap on proc "<<proc<<" :\n"<<endl;
  //      for (curr=slavegeom.begin(); curr!=slavegeom.end(); ++curr)
  //      {
  //        std::pair <std::multimap<int,double>::iterator, std::multimap<int,double>::iterator>
  //        range; range = slaveinterfaceelementidentificationmap.equal_range(curr->second->Id());
  //          std::cout<<curr->second->Id()<<" :";
  //          for (std::multimap<int,double>::iterator it=range.first; it!=range.second; ++it)
  //          {
  //            std::cout << ' ' << it->second;
  //          }
  //          std::cout<<"\n"<<endl;
  //      }
  //    }
  //    else
  //      mastercomm.Barrier();
  //  }

  //  int count=0;
  //  for (curr=slavegeom.begin(); curr!=slavegeom.end(); ++curr)
  //  {
  //    count++;
  //  }
  //  std::cout<<"count = "<< count << " "<<"on proc"<<mastercomm.MyPID()<<endl;
  ///////////////////////////////////////////////////////////////////////////////////
  ////////////                                                           ////////////
  ////////////                         MASTER                            ////////////
  ////////////                          Bulk                             ////////////
  ///////////////////////////////////////////////////////////////////////////////////

  // loop over procs
  for (int proc = 0; proc < mastercomm.NumProc(); proc++)
  {
    curr = mastergeom.begin();

    int parenteleid = -1;
    int parenteleowner = -1;
    int match = 0;
    int mastergeomsize = mastergeom.size();
    std::vector<int> sizelist(
        mastercomm.NumProc());  // how many master interface elements has each processor
    mastercomm.GatherAll(&mastergeomsize, sizelist.data(), 1);
    mastercomm.Barrier();  // wait for procs

    bool done;
    if (sizelist[proc])
      done = false;
    else
      done = true;  // no master interface elements on proc

    int counter = 0;
    while (not done)  // do until every master element on proc has been matched
    {
      match = 0;
      std::vector<double> masterloc;
      masterloc.assign(3, 0.0);

      if (masterdis.HaveDofs() == false)
      {
        FOUR_C_THROW("No degrees of freedom have been assigned to discretization, DaFUQ?!?!?");
      }

      // do for every master node
      if (proc == mastercomm.MyPID())
      {
        const int numnode = curr->second->num_node();

        for (int nodenum = 0; nodenum < numnode; nodenum++)
        {
          const Core::Nodes::Node* const* currmasternode = curr->second->Nodes();

          std::vector<double> temploc;
          temploc.assign(3, 0.0);
          temploc[0] = currmasternode[nodenum]->X()[0];
          temploc[1] = currmasternode[nodenum]->X()[1];
          temploc[2] = currmasternode[nodenum]->X()[2];
          for (int dim = 0; dim < 3; dim++)
          {
            masterloc[dim] = masterloc[dim] + temploc[dim];
          }

        }  // for every master node

        Teuchos::RCP<Core::Elements::FaceElement> bele =
            Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(curr->second);
        parenteleid = bele->parent_element()->Id();
        if (parenteleid == -1) FOUR_C_THROW("Couldn't get master parent element Id() ...");
        parenteleowner = bele->parent_element()->Owner();
        if (parenteleowner == -1) FOUR_C_THROW("Couldn't get master parent element Owner() ...");
      }

      mastercomm.Broadcast(&parenteleid, 1, proc);
      mastercomm.Broadcast(&parenteleowner, 1, proc);
      mastercomm.Broadcast(masterloc.data(), masterloc.size(), proc);

      mastercomm.Barrier();
      // match current master element
      // compare position to every element on interface slave side, every processor compares
      // masterloc of current master element of processor[proc]
      Teuchos::RCP<Core::Elements::Element> matchcurr = Teuchos::null;

      for (std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator scurr = slavegeom.begin();
           scurr != slavegeom.end(); ++scurr)
      {
        std::pair<std::multimap<int, double>::iterator, std::multimap<int, double>::iterator> range;
        range = slaveinterfaceelementidentificationmap.equal_range(scurr->second->Id());

        int dim = 0;
        match = 0;
        for (std::multimap<int, double>::iterator it = range.first; it != range.second; ++it)
        {
          double slaveloccomponent = it->second;
          double masterloccomponent = masterloc[dim];
          if (abs(slaveloccomponent - masterloccomponent) < 10e-3)
          {
            match += 1;
          }
          else
          {
            match += 0;
          }
          dim++;
        }

        if (match == 3)
        {
          matchcurr = scurr->second;
          break;
        }
      }  // end matching

      mastercomm.Barrier();

      for (int p = 0; p < mastercomm.NumProc(); ++p)
      {
        mastercomm.Barrier();
        if (mastercomm.MyPID() == p)
        {
          if (matchcurr != Teuchos::null)
          {
            // slave interface element is unique => key
            interfacefacingelementmap.insert(std::pair<int, int>(matchcurr->Id(), parenteleid));
          }
        }
      }

      mastercomm.Barrier();

      if (counter < (sizelist[mastercomm.MyPID()] - 1)) curr++;  // increment iterator

      if (counter != sizelist[proc]) counter++;
      if (counter == sizelist[proc])
        done = true;  // true on every proc
      else if (counter > sizelist[proc])
        FOUR_C_THROW("tried to match more master elements as are on proc ... ");

    }  // while iterating master elements on proc
    mastercomm.Barrier();
  }  // loop over procs

  int mymatchedelements = interfacefacingelementmap.size();
  int globalmatchedelements = 0;

  mastercomm.SumAll(&mymatchedelements, &globalmatchedelements, 1);

  if (mastercomm.MyPID() == 0)
    std::cout << "Could match " << globalmatchedelements << " " << slavedis.Name()
              << " interface elements to " << masterdis.Name() << " bulk elements." << std::endl;

  mastercomm.Barrier();  // wait for procs

  if (abs(globalslavegeomsize - globalmatchedelements) < 1e-6 and mastercomm.MyPID() == 0)
  {
    std::cout << "Setting up local interfacefacingelementmaps was successfull. \n" << std::endl;
  }
  else if (abs(globalslavegeomsize - globalmatchedelements) > 1e-3 and mastercomm.MyPID() == 0)
  {
    FOUR_C_THROW("ERROR: globalslavegeomsize != globalmatchedelements (%d,%d)", globalslavegeomsize,
        globalmatchedelements);
  }

  return;
}

/*---------------------------------------------------------------------------/
| Redistribute Interface (for parallel distr.)                        rauch  |
/---------------------------------------------------------------------------*/
void FPSI::Utils::redistribute_interface(Teuchos::RCP<Discret::Discretization> masterdis,
    Teuchos::RCP<const Discret::Discretization> slavedis, const std::string& condname,
    std::map<int, int>& interfacefacingelementmap)
{
  int printid = -1;

  std::map<int, int>::iterator mapcurr;
  std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator slaveelecurr;
  std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator masterelecurr;
  Core::Elements::Element* masterele = nullptr;

  Global::Problem* problem = Global::Problem::Instance();
  const Epetra_Comm& comm = problem->GetDis(masterdis->Name())->Comm();
  Teuchos::RCP<Epetra_Comm> rcpcomm = Teuchos::rcp(comm.Clone());

  int mymapsize = interfacefacingelementmap.size();
  int globalmapsize;
  std::vector<int> mapsizearray(comm.NumProc());
  comm.GatherAll(&mymapsize, mapsizearray.data(), 1);
  comm.SumAll(&mymapsize, &globalmapsize, 1);

  int counter = 0;
  for (int proc = 0; proc < comm.NumProc(); ++proc)
  {
    mapcurr = interfacefacingelementmap.begin();

    int done = 0;

    if (mapsizearray[proc] > 0)
      done = 0;
    else if (mapsizearray[proc] == 0)
      done = 1;
    else
      FOUR_C_THROW("weird size of processor local interfacefacingelementmap ... ");

    int HasMasterEle = 0;

    while (done == 0)
    {
      int slaveeleid = -1;
      int mastereleowner = -1;
      int mastereleid = -1;
      std::vector<int> slaveeleowners(comm.NumProc());
      std::vector<int> mastereleowners(comm.NumProc());

      // get master id
      if (comm.MyPID() == proc)
      {
        mastereleid = mapcurr->second;
        slaveeleid = mapcurr->first;
      }
      comm.Broadcast(&mastereleid, 1, proc);
      comm.Broadcast(&slaveeleid, 1, proc);

      if (masterdis->HaveGlobalElement(mastereleid))
      {
        masterele = masterdis->gElement(mastereleid);
        mastereleowner = masterele->Owner();

        if (masterele->Owner() != comm.MyPID())
        {
          masterele = nullptr;
          mastereleowner = -1;
        }
      }  // only the owner of the masterele has a pointer != nullptr and mastereleowner != -1

      comm.GatherAll(&mastereleowner, mastereleowners.data(), 1);

      for (int i = 0; i < comm.NumProc(); i++)
      {
        if (mastereleowners[i] != -1) mastereleowner = i;
      }  // now every processor knows the mastereleowner

      std::vector<int> procHasMasterEle(comm.NumProc());
      HasMasterEle = masterdis->HaveGlobalElement(mastereleid);
      comm.GatherAll(&HasMasterEle, procHasMasterEle.data(), 1);

      // ghost parent master element on master discretization of proc owning the matching slave
      // interface element
      const Epetra_Map colcopy = *(masterdis->ElementColMap());
      int myglobalelementsize = colcopy.NumMyElements();
      std::vector<int> myglobalelements(myglobalelementsize);
      colcopy.MyGlobalElements(myglobalelements.data());

      if (comm.MyPID() == proc and
          mastereleowner != proc)  // ghost master ele on owner of slave ele, but only if this proc
                                   // doesn't already own the masterele
      {
        if (colcopy.LID(mastereleid) == -1)  // if element not already in ElementColMap()
        {
          myglobalelements.push_back(mastereleid);
          myglobalelementsize = myglobalelementsize + 1;
        }
      }

      int globalsize;
      comm.SumAll(&myglobalelementsize, &globalsize, 1);
      Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(
          new Epetra_Map(globalsize, myglobalelementsize, myglobalelements.data(), 0, comm));

      if (mastereleid == printid)
      {
        std::cout << "mastereleowner" << mastereleowner << std::endl;
        std::cout << "slaveeleid" << slaveeleid << std::endl;
        std::cout << "counter" << counter << std::endl;
      }

      // Do the Redistribution
      int before = 0;
      int after = 0;
      if (procHasMasterEle[proc] == 0)
      {
        if (comm.MyPID() == proc)
        {
          // std::cout<<counter<<" --Before: Have GID "<<mastereleid<<" =
          // "<<masterdis->HaveGlobalElement(mastereleid)<<" on proc "<<slaveeleowner<<endl;
          before = masterdis->HaveGlobalElement(mastereleid);
        }
        comm.Barrier();
        masterdis->ExtendedGhosting(*newelecolmap, true, false, true, true);
        if (comm.MyPID() == proc)
        {
          // std::cout<<counter<<" --After: Have GID "<<mastereleid<<" =
          // "<<masterdis->HaveGlobalElement(mastereleid)<<" on proc "<<slaveeleowner<<endl;
          after = masterdis->HaveGlobalElement(mastereleid);
          if (after == 0 and before == 0)
            FOUR_C_THROW("Element with gid=%d has not been redistributed ! ", mastereleid);
        }
      }

      if (comm.MyPID() == proc and mapcurr != interfacefacingelementmap.end()) ++mapcurr;
      if ((comm.MyPID() == proc and mapcurr == interfacefacingelementmap.end())) done = 1;

      if (counter == globalmapsize) done = 1;
      comm.Broadcast(&done, 1, proc);

      if (done == 0) counter++;

    }  // for all elements of interfacefacingelementmap
    comm.Barrier();
  }  // for all procs

  if (comm.MyPID() == 0)
    std::cout << "\n EXTENDEDGHOSTING: checked " << counter
              << " elements in interfacefacingelementmap ... \n"
              << std::endl;
  return;
}

/*---------------------------------------------------------------------------/
| Setup Interface Map (for parallel distr.)                           rauch  |
/---------------------------------------------------------------------------*/
void FPSI::Utils::SetupInterfaceMap(const Epetra_Comm& comm,
    Teuchos::RCP<Discret::Discretization> structdis,
    Teuchos::RCP<Discret::Discretization> porofluiddis,
    Teuchos::RCP<Discret::Discretization> fluiddis, Teuchos::RCP<Discret::Discretization> aledis)
{
  poro_fluid_fluid_interface_map_ = Teuchos::rcp(new std::map<int, int>);
  fluid_poro_fluid_interface_map_ = Teuchos::rcp(new std::map<int, int>);

  setup_local_interface_facing_element_map(
      *fluiddis, *porofluiddis, "fpsi_coupling", *poro_fluid_fluid_interface_map_);
  setup_local_interface_facing_element_map(
      *porofluiddis, *fluiddis, "fpsi_coupling", *fluid_poro_fluid_interface_map_);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::UTILS::MapExtractor::Setup(
    const Discret::Discretization& dis, bool withpressure, bool overlapping)
{
  const int ndim = Global::Problem::Instance()->NDim();
  Core::Conditions::MultiConditionSelector mcs;
  mcs.SetOverlapping(overlapping);  // defines if maps can overlap
  mcs.AddSelector(Teuchos::rcp(
      new Core::Conditions::NDimConditionSelector(dis, "FSICoupling", 0, ndim + withpressure)));
  mcs.AddSelector(Teuchos::rcp(
      new Core::Conditions::NDimConditionSelector(dis, "fpsi_coupling", 0, ndim + withpressure)));
  mcs.SetupExtractor(dis, *dis.dof_row_map(), *this);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::UTILS::MapExtractor::Setup(
    Teuchos::RCP<const Epetra_Map>& additionalothermap, const FPSI::UTILS::MapExtractor& extractor)
{
  // build the new othermap
  std::vector<Teuchos::RCP<const Epetra_Map>> othermaps;
  othermaps.push_back(additionalothermap);
  othermaps.push_back(extractor.OtherMap());

  if (Core::LinAlg::MultiMapExtractor::IntersectMaps(othermaps)->NumGlobalElements() > 0)
    FOUR_C_THROW("Failed to add dofmap of foreign discretization to OtherMap. Detected overlap.");

  Teuchos::RCP<const Epetra_Map> mergedothermap =
      Core::LinAlg::MultiMapExtractor::MergeMaps(othermaps);

  // the vector of maps for the new map extractor consists of othermap at position 0
  // followed by the maps of conditioned DOF
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  // append the merged other map at first position
  maps.push_back(mergedothermap);

  // append the condition maps subsequently
  for (int i = 1; i < extractor.NumMaps(); ++i) maps.push_back(extractor.Map(i));

  // merge
  Teuchos::RCP<const Epetra_Map> fullmap = Core::LinAlg::MultiMapExtractor::MergeMaps(maps);

  Core::LinAlg::MultiMapExtractor::Setup(*fullmap, maps);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int>> FPSI::UTILS::MapExtractor::conditioned_element_map(
    const Discret::Discretization& dis) const
{
  Teuchos::RCP<std::set<int>> condelements =
      Core::Conditions::conditioned_element_map(dis, "FSICoupling");
  Teuchos::RCP<std::set<int>> condelements2 =
      Core::Conditions::conditioned_element_map(dis, "fpsi_coupling");
  std::copy(condelements2->begin(), condelements2->end(),
      std::inserter(*condelements, condelements->begin()));
  return condelements;
}

FOUR_C_NAMESPACE_CLOSE
