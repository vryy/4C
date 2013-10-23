/*----------------------------------------------------------------------*/
/*!
 \file fpsi_utils.cpp

 \brief  Utility Methods For Fluid Porous Structure Interaction Problems
 <pre>
   Maintainer: Rauch Andreas
               rauch@lnm.mw.tum.de
 </pre>
 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  rauch 11/12 |
 *----------------------------------------------------------------------*/
//GENERAL
#include <Teuchos_TimeMonitor.hpp>
#include <sstream>

//FPSI includes
#include "fpsi.H"
#include "fpsi_dyn.H"
#include "fpsi_utils.H"
#include "fpsi_monolithic.H"
#include "fpsi_partitioned.H"
#include "fpsi_robinneumann.H"
#include "fpsi_monolithic_plain.H"

//POROELAST includes
#include "../drt_poroelast/poro_base.H"
#include "../drt_poroelast/poroelast_utils.H"

//FSI includes
#include "../drt_fsi/fsi_utils.H"

#include "../drt_fsi/fsi_algorithm.H"
#include "../drt_ale/ale_utils_clonestrategy.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"

//INPAR includes
#include "../drt_inpar/inpar_fpsi.H"
#include "../drt_inpar/inpar_solver.H"

//lib includes
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parmetis.H"

//ADAPTER includes
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

//LINALG includes
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"

//IO includes
#include "../drt_io/io_control.H"

//FLUID includes
#include "../drt_fluid_ele/fluid_ele.H"

//STRUCTURE includes
#include "../drt_structure/stru_aux.H"


Teuchos::RCP<FPSI::UTILS> FPSI::UTILS::instance_;

Teuchos::RCP<FPSI::UTILS> FPSI::UTILS::Instance()
{
  if(instance_==Teuchos::null)
    instance_ = Teuchos::rcp(new UTILS());

  return instance_;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    Setup Discretization    >>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Teuchos::RCP<FPSI::FPSI_Base > FPSI::UTILS::SetupDiscretizations(const Epetra_Comm& comm,
                                                                 const Teuchos::ParameterList& fpsidynparams,
                                                                 const Teuchos::ParameterList& poroelastdynparams)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  Fluid_PoroFluid_InterfaceMap = Teuchos::rcp(new std::map<int,int>);
  PoroFluid_Fluid_InterfaceMap = Teuchos::rcp(new std::map<int,int>);

  //1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis     = problem->GetDis("structure");
  if(structdis == Teuchos::null)
    dserror(" !!! structdis empty !!! Awwww MAAAAN !!!");
  Teuchos::RCP<DRT::Discretization> porofluiddis  = problem->GetDis("porofluid");
  if(porofluiddis == Teuchos::null)
    dserror(" !!! porofluiddis empty !!! Awwww MAAAAN !!!");
  Teuchos::RCP<DRT::Discretization> fluiddis      = problem->GetDis("fluid");
  if(fluiddis== Teuchos::null)
    dserror(" !!! fluiddis empty !!! Awwww MAAAAN !!!");
  Teuchos::RCP<DRT::Discretization> aledis        = problem->GetDis("ale");
  if(aledis == Teuchos::null)
    dserror(" !!! aledis empty !!! Awwww MAAAAN !!!");

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
 diagonal blocks (from northwest to southeast): porostructure - porofluid - fee fluid - ale.

Thats nice....

This would be corrupted however if we were doing the following:

structdis    -> FillComplete()
porofluiddis -> FillComplete()
fluiddis     -> FillComplete()
aledis       -> FillComplete()

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
  structdis->FillComplete(true,true,true);
  if (structdis->NumGlobalNodes() == 0)
  {
     dserror("Structure discretization is empty! Damn it!");
  }

  porofluiddis->FillComplete(true,true,true);

  if (porofluiddis->NumGlobalNodes() == 0)
  {
    // Create the fluid discretization within the porous medium
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroelastCloneStrategy>(structdis,porofluiddis);
  }
  else
  {
    dserror("Porofluiddis already exists! That's quite impossible!");
  }

  fluiddis->FillComplete(true,true,true);

  aledis->FillComplete(true,true,true);

  //3.- Create ALE elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0) // ALE discretization still empty
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
  }
  else  // ALE discretization already filled
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis,aledis))
      dserror("Fluid and ALE nodes have the same node numbers. "
              "This it not allowed since it causes problems with Dirichlet BCs. "
              "Use the ALE cloning functionality or ensure non-overlapping node numbering!"
             );
  }

  PoroFluid_Fluid_InterfaceMap = Teuchos::rcp(new std::map<int,int> );
  Fluid_PoroFluid_InterfaceMap = Teuchos::rcp(new std::map<int,int> );

  SetupLocalInterfaceFacingElementMap(*fluiddis,*porofluiddis,"FSICoupling",*PoroFluid_Fluid_InterfaceMap);
  SetupLocalInterfaceFacingElementMap(*porofluiddis,*fluiddis,"FSICoupling",*Fluid_PoroFluid_InterfaceMap);
  if(comm.NumProc() > 1 and problem->Restart() == 0)
  {
    RedistributeInterface(fluiddis,*porofluiddis,"FSICoupling",*PoroFluid_Fluid_InterfaceMap);
    RedistributeInterface(aledis,*porofluiddis,"FSICoupling",*PoroFluid_Fluid_InterfaceMap);
    RedistributeInterface(porofluiddis,*fluiddis,"FSICoupling",*Fluid_PoroFluid_InterfaceMap);
    RedistributeInterface(structdis,*fluiddis,"FSICoupling",*Fluid_PoroFluid_InterfaceMap);
  }

    //4.- get coupling algorithm
    Teuchos::RCP<FPSI::FPSI_Base> fpsi_algo = Teuchos::null;
    int coupling = DRT::INPUT::IntegralValue<int>(fpsidynparams,"COUPALGO");
    switch (coupling)
      {
        case fpsi_monolithic_plain:
        {
          fpsi_algo = Teuchos::rcp(new FPSI::Monolithic_Plain(comm,fpsidynparams,poroelastdynparams));
          break;
        } // case monolithic
        case partitioned:
        {
          dserror("Partitioned solution scheme not implemented for FPSI, yet. "
                  "Make sure that the parameter COUPALGO is set to 'fpsi_monolithic_plain', "
                  "and the parameter PARITIONED is set to 'monolithic'. ");
          INPAR::FPSI::PartitionedCouplingMethod method;
          method = DRT::INPUT::IntegralValue<INPAR::FPSI::PartitionedCouplingMethod>(fpsidynparams,"PARTITIONED");
          if(method == INPAR::FPSI::RobinNeumann )
          {
            fpsi_algo = Teuchos::rcp(new FPSI::RobinNeumann(comm,fpsidynparams));
          }
          else
          {
            dserror("Only RobinNeumann algorithm available for partitioned FPSI !!!\n"
                    "Set 'PARTITIONED' to 'RobinNeumann' in input file.");
          }
          break;
        }
      } // switch(coupling)

    return fpsi_algo;
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<          FINISHED          >>>>>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    Setup Discretization    >>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void FPSI::UTILS::SetupLocalInterfaceFacingElementMap(DRT::Discretization& masterdis, const DRT::Discretization& slavedis, const std::string& condname, std::map<int,int>& interfacefacingelementmap)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  const Epetra_Comm& mastercomm = problem->GetDis(masterdis.Name())->Comm();


  DRT::Condition* slavecond = slavedis.GetCondition(condname);
  std::map<int,RCP<DRT::Element> >& slavegeom = slavecond->Geometry();

  DRT::Condition* mastercond = masterdis.GetCondition(condname);
  std::map<int,RCP<DRT::Element> >& mastergeom = mastercond->Geometry();

  std::map<int,RCP<DRT::Element> >::iterator curr;

  std::multimap<int,double> slaveinterfaceelementidentificationmap;


  ///////////////////////////////////////////////////////////////////////////////////
  ////////////                                                           ////////////
  ////////////                         SLAVE                             ////////////
  ////////////                       Interface                           ////////////
  ///////////////////////////////////////////////////////////////////////////////////
  int slavegeomsize = slavegeom.size();
  int globalslavegeomsize;
  mastercomm.SumAll(&slavegeomsize,&globalslavegeomsize,1);

  // do for every slave interface element
  for (curr=slavegeom.begin(); curr!=slavegeom.end(); ++curr)
   {
    //std::cout<<"Proc "<<mastercomm.MyPID()<<" owns slavegeom "<<curr->first<<endl;

     std::vector<double> slaveloc;
     slaveloc.assign(3,0.0);

//       std::cout<<"Current Slave Interface Element ID: "<<curr->second->Id()<<endl;

     if(slavedis.HaveDofs() == false)
     {
       dserror("No degrees of freedom have been assigned to discretization");
     }

     // do for every interface slave node
     for(int nodenum=0; nodenum < curr->second->NumNode(); nodenum++)
     {
       const DRT::Node*const* currslavenode  = curr->second->Nodes();

       std::vector<double> temploc;
       temploc.assign(3,0.0);
       temploc[0] = currslavenode[nodenum]->X()[0];
       temploc[1] = currslavenode[nodenum]->X()[1];
       temploc[2] = currslavenode[nodenum]->X()[2];
       for(int dim=0; dim<3; dim++)
       {
         slaveloc[dim]  = slaveloc[dim] + temploc[dim];
       }
     }// for every slave node

     slaveinterfaceelementidentificationmap.insert(std::pair<int,double>(curr->second->Id(),slaveloc[0]));
     slaveinterfaceelementidentificationmap.insert(std::pair<int,double>(curr->second->Id(),slaveloc[1]));
     slaveinterfaceelementidentificationmap.insert(std::pair<int,double>(curr->second->Id(),slaveloc[2]));

   } // do for every global slave element

  mastercomm.Barrier(); // wait for procs

//  // printout
//  for(int proc = 0; proc < (mastercomm.NumProc()+1); proc++)
//  {
//    if(mastercomm.MyPID() == proc)
//    {
//      std::cout<<"\n slaveinterfaceelementidentificationmap on proc "<<proc<<" :\n"<<endl;
//      for (curr=slavegeom.begin(); curr!=slavegeom.end(); ++curr)
//      {
//        std::pair <std::multimap<int,double>::iterator, std::multimap<int,double>::iterator> range;
//        range = slaveinterfaceelementidentificationmap.equal_range(curr->second->Id());
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
  for(int proc = 0; proc < mastercomm.NumProc(); proc++)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidBoundary> fbele = Teuchos::null;
    curr=mastergeom.begin();
//    std::cout<<curr->second->Id()<<" proc :" <<mastercomm.MyPID()<<endl;
    int parenteleid=-1;
    int parenteleowner=-1;
    int match = 0;
    int mastergeomsize = mastergeom.size();
    int sizelist[mastercomm.NumProc()]; // how many master interface elements has each processor
    mastercomm.GatherAll(&mastergeomsize,&sizelist[0],1);
    mastercomm.Barrier(); // wait for procs
//    if(mastercomm.MyPID() == proc)
//    {
//      for(int i=0; i<mastercomm.NumProc();i++)
//        std::cout<<sizelist[i]<<endl;
//    }

    bool done;
    if(sizelist[proc])
      done = false;
    else
      done = true; // no master interface elements on proc

    int counter = 0;
    while(not done) // do until every master element on proc has been matched
    {
      match = 0;
      std::vector<double> masterloc;
      masterloc.assign(3,0.0);

      if(masterdis.HaveDofs() == false)
      {
        dserror("No degrees of freedom have been assigned to discretization, DaFUQ?!?!?");
      }

      // do for every master node
      if(proc == mastercomm.MyPID())
      {
        int numnode = curr->second->NumNode();

        if(numnode != 4)
        {
          std::cout<<"numnode = "<<numnode<<std::endl;
          std::cout<<"Proc :    "<< mastercomm.MyPID()<<std::endl;
          std::cout<<"variable proc    : "<< proc <<std::endl;
          std::cout<<"variable counter : "<< counter <<std::endl;
          dserror("Only Hex8 parent elements supported so far. Interface geometry elements must have 4 nodes.");
        }
        for(int nodenum=0; nodenum < numnode; nodenum++)
        {
          const DRT::Node*const* currmasternode  = curr->second->Nodes();

          std::vector<double> temploc;
          temploc.assign(3,0.0);
          temploc[0] = currmasternode[nodenum]->X()[0];
          temploc[1] = currmasternode[nodenum]->X()[1];
          temploc[2] = currmasternode[nodenum]->X()[2];
          for(int dim=0; dim<3; dim++)
          {
            masterloc[dim]  = masterloc[dim] + temploc[dim];
          }

        } // for every master node

        Teuchos::RCP<DRT::ELEMENTS::FluidBoundary> fbele = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidBoundary>(curr->second);
        parenteleid = fbele->ParentElement()->Id();
        if(parenteleid == -1)
          dserror("Couldn't get master parent element Id() ...");
        parenteleowner = fbele->ParentElement()->Owner();
        if(parenteleid == -1)
          dserror("Couldn't get master parent element Owner() ...");

      }

      mastercomm.Broadcast(&parenteleid,1,proc);
      mastercomm.Broadcast(&parenteleowner,1,proc);
      mastercomm.Broadcast(&masterloc[0],masterloc.size(),proc);


//      for(int j=0; j<3; j++)
//      {
//        std::cout<<" "<<masterloc[j]<<" "<<"on proc "<<mastercomm.MyPID()<<endl;;
//      }

      mastercomm.Barrier();

      // match current master element
      // compare position to every element on interface slave side, every processor compares masterloc of current master element of processor[proc]
      Teuchos::RCP<DRT::Element> matchcurr = Teuchos::null;

      for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator scurr=slavegeom.begin(); scurr!=slavegeom.end(); ++scurr)
      {
        std::pair <std::multimap<int,double>::iterator, std::multimap<int,double>::iterator> range;
        range = slaveinterfaceelementidentificationmap.equal_range(scurr->second->Id());

        int dim = 0;
        match = 0;
        for (std::multimap<int,double>::iterator it=range.first; it!=range.second; ++it)
        {
          double slaveloccomponent  = it->second;
          double masterloccomponent = masterloc[dim];
          if(abs(slaveloccomponent - masterloccomponent) < 10e-3)
          {
            match += 1;
          }
          else
          {
            match += 0;
          }
          dim++;
        }

        if(match == problem->NDim())
        {
          matchcurr = scurr->second;
          break;
        }
      } // end matching

      mastercomm.Barrier();

      if(matchcurr != Teuchos::null)
      {
        //std::cout<<"FOUND MATCH ON PROC "<<mastercomm.MyPID()<<" for slaveelement "<<matchcurr->Id()<<" and master element "<<parenteleid<<" !!! ... "<<endl;
        interfacefacingelementmap.insert(std::pair<int,int>(matchcurr->Id(),parenteleid)); // slave interface element is unique => key
      }

      mastercomm.Barrier();

      if(counter<(sizelist[mastercomm.MyPID()]-1))
        curr++; // increment iterator

      if(counter != sizelist[proc])
        counter++;
      if(counter == sizelist[proc])
        done = true; // true on every proc
      else if(counter > sizelist[proc])
        dserror("tried to match more master elements as are on proc ... ");


    } // while iterating master elements on proc
    mastercomm.Barrier();
  } // loop over procs

 int mymatchedelements = interfacefacingelementmap.size();
 int globalmatchedelements = 0;

 mastercomm.SumAll(&mymatchedelements,&globalmatchedelements,1);

 if(mastercomm.MyPID()==0)
   std::cout<<"Could match "<<globalmatchedelements<<" "<<slavedis.Name()<<" interface elements to "<<masterdis.Name()<<" bulk elements."<<std::endl;

// std::cout << masterdis.Name() << " => " << slavedis.Name() <<endl;
// for ( std::map<int,int>::iterator it=interfacefacingelementmap.begin(); it!=interfacefacingelementmap.end(); ++it)
//     std::cout << it->first << " => " << it->second <<endl;

 mastercomm.Barrier(); // wait for procs

 if(abs(globalslavegeomsize - globalmatchedelements) < 1e-6 and mastercomm.MyPID() == 0)
   std::cout<<"Setting up local interfacefacingelementmaps was successfull. \n"<<std::endl;
 else if(abs(globalslavegeomsize - globalmatchedelements) > 1e-3 and mastercomm.MyPID() == 0)
   dserror("ERROR");
return;
}

void FPSI::UTILS::RedistributeInterface(Teuchos::RCP<DRT::Discretization> masterdis, const DRT::Discretization& slavedis, const std::string& condname, std::map<int,int>& interfacefacingelementmap)
{

  int printid = -1;

  std::map<int,int>::iterator mapcurr;
  std::map<int,RCP<DRT::Element> >::iterator slaveelecurr;
  std::map<int,RCP<DRT::Element> >::iterator masterelecurr;
  DRT::Element* masterele = NULL;

  DRT::Problem* problem = DRT::Problem::Instance();
  const Epetra_Comm& comm = problem->GetDis(masterdis->Name())->Comm();
  Teuchos::RCP<Epetra_Comm> rcpcomm = Teuchos::rcp(comm.Clone());

  //DRT::Condition* slavecond = slavedis.GetCondition(condname);
  //std::map<int,RCP<DRT::Element> >& slavegeom = slavecond->Geometry();

  int mymapsize = interfacefacingelementmap.size();
  int globalmapsize;
  std::vector<int> mapsizearray(comm.NumProc());
  comm.GatherAll(&mymapsize,&mapsizearray[0],1);
  comm.SumAll(&mymapsize,&globalmapsize,1);
  //std::cout<<"Proc "<<comm.MyPID()<<" : "<<mymapsize<<endl;

int counter = 0;
for(int proc=0; proc<comm.NumProc(); ++proc)
{
  mapcurr=interfacefacingelementmap.begin();

  int done;

  if(mapsizearray[proc] > 0)
    done = 0;
  else if(mapsizearray[proc] == 0)
    done = 1;
  else
    dserror("weird size of processor local interfacefacingelementmap ... ");

  int HasMasterEle = 0;

  while(done == 0)
  {
    int slaveeleid=-1;
    int mastereleowner=-1;
    int mastereleid=-1;
    std::vector<int> slaveeleowners (comm.NumProc());
    std::vector<int> mastereleowners(comm.NumProc());

    // get master id
    if(comm.MyPID()==proc)
    {
      mastereleid = mapcurr->second;
      slaveeleid  = mapcurr->first;
    }
    comm.Broadcast(&mastereleid,1,proc);
    comm.Broadcast(&slaveeleid,1,proc);

    if(masterdis->HaveGlobalElement(mastereleid))
    {
      masterele = masterdis->gElement(mastereleid);
      mastereleowner = masterele->Owner();

      if(masterele->Owner() != comm.MyPID())
      {
        masterele = NULL;
        mastereleowner = -1;
      }
    } // only the owner of the masterele has a pointer != NULL and mastereleowner != -1

    comm.GatherAll(&mastereleowner,&mastereleowners[0],1);

    for(int i=0;i<comm.NumProc();i++)
    {
      if(mastereleowners[i] != -1)
        mastereleowner = i;
    } // now every processor knows the mastereleowner


    std::vector<int> procHasMasterEle(comm.NumProc());
    HasMasterEle = masterdis->HaveGlobalElement(mastereleid);
    comm.GatherAll(&HasMasterEle,&procHasMasterEle[0],1);


    // ghost parent master element on master discretization of proc owning the matching slave interface element
    const Epetra_Map colcopy = *(masterdis->ElementColMap());
    int myglobalelementsize = colcopy.NumMyElements();
    std::vector<int> myglobalelements(myglobalelementsize);
    colcopy.MyGlobalElements(&myglobalelements[0]);

    if(comm.MyPID() == proc and mastereleowner != proc) // ghost master ele on owner of slave ele, but only if this proc doesn't already own the masterele
    {
      if(colcopy.LID(mastereleid) == -1) // if element not already in ElementColMap()
      {
        myglobalelements.push_back(mastereleid);
        myglobalelementsize = myglobalelementsize + 1;
      }
    }

    int globalsize;
    comm.SumAll(&myglobalelementsize,&globalsize,1);
    Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(new Epetra_Map(globalsize,myglobalelementsize,&myglobalelements[0],0,comm));


    if(mastereleid == printid)
    {
      std::cout<<"mastereleowner"<<mastereleowner<<std::endl;
      std::cout<<"slaveeleid"<<slaveeleid<<std::endl;
      std::cout<<"counter"<<counter<<std::endl;
    }

    // Do the Redistribution
    int before=0;
    int after=0;
    if(procHasMasterEle[proc] == 0)
    {
      if(comm.MyPID() == proc)
      {
        //std::cout<<counter<<" --Before: Have GID "<<mastereleid<<" = "<<masterdis->HaveGlobalElement(mastereleid)<<" on proc "<<slaveeleowner<<endl;
        before = masterdis->HaveGlobalElement(mastereleid);
      }
      comm.Barrier();
      masterdis->ExtendedGhosting(*newelecolmap,true,false,true,true);
      if(comm.MyPID() == proc)
      {
        //std::cout<<counter<<" --After: Have GID "<<mastereleid<<" = "<<masterdis->HaveGlobalElement(mastereleid)<<" on proc "<<slaveeleowner<<endl;
        after = masterdis->HaveGlobalElement(mastereleid);
        if(after == 0 and before == 0)
          dserror("Element with gid=%d has not been redistributed ! ",mastereleid);
      }
    }

    if(comm.MyPID()==proc and mapcurr != interfacefacingelementmap.end())
      ++mapcurr;
    else if((comm.MyPID()==proc and mapcurr == interfacefacingelementmap.end()))
      done = 1;


    if(counter==globalmapsize)
      done=1;
    comm.Broadcast(&done,1,proc);

    if(done==0)
      counter++;

  } // for all elements of interfacefacingelementmap
  comm.Barrier();
} // for all procs

if(comm.MyPID()==0)
  std::cout<<"\n EXTENDEDGHOSTING: checked "<<counter<<" elements in interfacefacingelementmap ... \n"<<std::endl;
  return;
}

