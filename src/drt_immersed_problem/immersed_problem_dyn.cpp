/*!----------------------------------------------------------------------
\file immersed_problem_dyn.cpp

\brief global algorithm control class for all immersed and cell migration algorithms

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

\level 2

*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "immersed_problem_dyn.H"
#include "immersed_base.H"
#include "immersed_node.H"
#include "immersed_partitioned.H"
#include "immersed_partitioned_fsi.H"
#include "immersed_partitioned_confine_cell.H"
#include "immersed_partitioned_cellmigration.H"
#include "ssi_partitioned_2wc_biochemomechano.H"
#include "immersed_partitioned_adhesion_traction.H"
#include "ssi_partitioned_2wc_protrusionformation.H"
#include "immersed_partitioned_fsi_dirichletneumann.H"
#include "immersed_partitioned_protrusion_formation.H"
#include "immersed_partitioned_flow_cell_interaction.H"
#include "immersed_partitioned_fsi_dirichletneumann_ale.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_inpar/inpar_immersed.H"
#include "../drt_inpar/inpar_cell.H"
#include "../drt_inpar/inpar_ssi.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_ale/ale_utils_clonestrategy.H"
#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_poroelast/poroelast_utils_setup.H"
#include "../drt_poroelast/poroelast_monolithic.H"
#include "../drt_poroelast/poro_scatra_part_2wc.H"
#include "../linalg/linalg_utils.H"
#include "../drt_scatra/scatra_dyn.H"

void immersed_problem_drt()
{
  // get pointer to global problem
  DRT::Problem* problem = DRT::Problem::Instance();
  // get communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();
  // get parameterlist for immersed method
  const Teuchos::ParameterList& immersedmethodparams = problem->ImmersedMethodParams();
  // choose algorithm
  int coupling = DRT::INPUT::IntegralValue<int>(immersedmethodparams,"COUPALGO");
  int scheme   = DRT::INPUT::IntegralValue<int>(immersedmethodparams,"SCHEME");
  // declare general ParameterList that can be handed into Algorithm
  Teuchos::ParameterList params;


  switch (coupling)
  {
  case INPAR::IMMERSED::partitioned:
  {
    switch(DRT::Problem::Instance()->ProblemType())
    {
    case prb_immersed_fsi:
    {
      // SAFETY FIRST
      {
        // check if INODE is defined in input file
        int gid = problem->GetDis("fluid")->ElementRowMap()->GID(0);
        IMMERSED::ImmersedNode* inode =
            dynamic_cast<IMMERSED::ImmersedNode* >((problem->GetDis("fluid")->gElement(gid)->Nodes()[0]));

            if(inode == NULL)
              dserror("dynamic cast from Node to ImmersedNode failed.\n"
                      "Make sure you defined INODE instead of NODE in your input file.");
      }

      {
        // check if structural predictor ConstDisVelAcc is chosen in input file
        if(problem->StructuralDynamicParams().get<std::string>("PREDICT") != "ConstDisVelAcc")
          dserror("Invalid structural predictor for immersed fsi!\n"
                  "Choose ConstDisVelAcc as predictor in ---STRUCTURAL DYNAMIC section.\n"
                  "Structural state projected onto fluid in new time step should be the same as in previous time step.");
      }

      // fill discretizations
      problem->GetDis("structure")->FillComplete();
      problem->GetDis("fluid")    ->FillComplete();

      Teuchos::RCP<IMMERSED::ImmersedPartitionedFSIDirichletNeumann> algo = Teuchos::null;
      if(scheme == INPAR::IMMERSED::dirichletneumann)
        algo = Teuchos::rcp(new IMMERSED::ImmersedPartitionedFSIDirichletNeumann(comm));
      else
      {
        algo = Teuchos::null;
        dserror("unknown coupling scheme");
      }

      // ghost structure redundantly on all procs
      algo->CreateGhosting(problem->GetDis("structure"));

      // setup algo
      algo->Setup();

      // PARTITIONED FSI ALGORITHM

      // read restart step
      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        algo->ReadRestart(restart);
      }

      // additional setup for structural search tree, etc.
      algo->SetupStructuralDiscretization();

      algo->Timeloop(algo);

      if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
      {
        Teuchos::TimeMonitor::summarize();
        Teuchos::TimeMonitor::zeroOutTimers();
      }

      // create result tests for single fields
      DRT::Problem::Instance()->AddFieldTest(algo->MBFluidField()->CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(algo->StructureField()->CreateFieldTest());

      // do the actual testing
      DRT::Problem::Instance()->TestAll(comm);

      break;
    }// case prb_immersed_fsi
    case prb_immersed_ale_fsi:
    {
      // SAFETY FIRST
      {
        // check if INODE is defined in input file
        int gid = problem->GetDis("fluid")->ElementRowMap()->GID(0);
        if(gid!=-1)
        {
          IMMERSED::ImmersedNode* inode =
              dynamic_cast<IMMERSED::ImmersedNode* >((problem->GetDis("fluid")->gElement(gid)->Nodes()[0]));

            if(inode == NULL)
              dserror("dynamic cast from Node to ImmersedNode failed.\n"
                      "Make sure you defined INODE instead of NODE in your input file.");
        }
      }

      {
        // check if structural predictor ConstDisVelAcc is chosen in input file
        if(problem->StructuralDynamicParams().get<std::string>("PREDICT") != "ConstDisVelAcc")
          dserror("Invalid structural predictor for immersed fsi!\n"
                  "Choose ConstVel as predictor in ---STRUCTURAL DYNAMIC section.\n"
                  "Structural state projected onto fluid in new time step should be the same as in previous time step.");
      }

      // fill discretizations
      problem->GetDis("structure")->FillComplete();
      problem->GetDis("fluid")    ->FillComplete();
      problem->GetDis("ale")      ->FillComplete();

      // get discretizations
      Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
      Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

      // create ale elements if the ale discretization is empty
      if (aledis->NumGlobalNodes()==0)
      {
        DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
        // setup material in every ALE element
        Teuchos::ParameterList params;
        params.set<std::string>("action", "setup_material");
        aledis->Evaluate(params);
      }
      else  // filled ale discretization
      {
        if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis,aledis))
          dserror("Fluid and ALE nodes have the same node numbers. "
                  "This it not allowed since it causes problems with Dirichlet BCs. "
                  "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
      }

      // create algorithm
      Teuchos::RCP<IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE> algo = Teuchos::null;
      if(scheme == INPAR::IMMERSED::dirichletneumann)
        algo = Teuchos::rcp(new IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE(comm));
      else
      {
        algo = Teuchos::null;
        dserror("unknown coupling scheme");
      }

      // setup algo
      algo->Setup();

      // additional setup of structural search tree, etc.
      algo->SetupStructuralDiscretization();

      // PARTITIONED FSI ALGORITHM
      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        algo->ReadRestart(restart);
      }

      algo->Timeloop(algo);

      if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
      {
        Teuchos::TimeMonitor::summarize();
        Teuchos::TimeMonitor::zeroOutTimers();
      }

      // create result tests for single fields
      DRT::Problem::Instance()->AddFieldTest(algo->MBFluidField()->CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(algo->StructureField()->CreateFieldTest());

      // do the actual testing
      DRT::Problem::Instance()->TestAll(comm);

      break;
    }// case prb_immersed_fsi
    case prb_immersed_cell:
    {
      CellMigrationControlAlgorithm();

      break;
    }// case prb_cell_migration
    default:
    {
      dserror("no valid problem type specified");
      break;
    }// default
    break;
    }//switch problemtype
    break;
  } // case partitioned (default)
  case INPAR::IMMERSED::monolithic:
  {
    dserror("Monolithic solution scheme not implemented for immersed problems, yet.\n "
            "Make sure that the parameter COUPALGO is set to 'partitioned'");
    break;
  }// case monolithic
  }// end switch(coupling)

}// immersed_problem_drt()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CellMigrationControlAlgorithm()
{

  // get pointer to global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // get parameterlist for immersed method
  const Teuchos::ParameterList& immersedmethodparams = problem->ImmersedMethodParams();

  // declare general ParameterList that can be handed into Algorithm
  Teuchos::ParameterList params;

  // get parameterlist for cell migration
  const Teuchos::ParameterList& cellmigrationparams = problem->CellMigrationParams();

  // extract the simulation type
  int simtype = DRT::INPUT::IntegralValue<int>(cellmigrationparams,"SIMTYPE");

  //use PoroelastImmersedCloneStrategy to build FluidPoroImmersed elements
  POROELAST::UTILS::SetupPoroScatraDiscretizations
  <POROELAST::UTILS::PoroelastImmersedCloneStrategy,POROELAST::UTILS::PoroScatraCloneStrategy>();

  {
    // check if INODE is defined in input file
    int gid = problem->GetDis("porofluid")->ElementRowMap()->GID(0);
    if(gid!=-1)
    {
      IMMERSED::ImmersedNode* inode =
          dynamic_cast<IMMERSED::ImmersedNode* >((problem->GetDis("porofluid")->gElement(gid))->Nodes()[0]);

      if(inode == NULL)
        dserror("dynamic cast from Node to ImmersedNode failed.\n"
            "Make sure you defined INODE instead of NODE in your input file.");
    }
  }

  // assign degrees of freedom, initialize elements and do boundary conditions on cell
  problem->GetDis("cell")->FillComplete(true,true,true);

  // pointer to field cell structure
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapperImmersed> cellstructure = Teuchos::null;

  // pointer to cell subproblem (structure-scatra interaction)
  Teuchos::RCP<SSI::SSI_Part2WC_PROTRUSIONFORMATION> cellscatra_subproblem = Teuchos::null;

  // check if cell is supposed to have intracellular biochchemical reaction capabilities
  bool ssi_cell = DRT::INPUT::IntegralValue<int>(problem->CellMigrationParams(),"SSI_CELL");
  if(ssi_cell and simtype==INPAR::CELL::sim_type_pureProtrusionFormation)
  {
    // get coupling algorithm from ---SSI CONTROL section
    const INPAR::SSI::SolutionSchemeOverFields coupling
    = DRT::INPUT::IntegralValue<INPAR::SSI::SolutionSchemeOverFields>(problem->SSIControlParams(),"COUPALGO");

    switch(coupling)
    {
    case INPAR::SSI::ssi_IterStagg: // create SSI_Part2WC_PROTRUSIONFORMATION subproblem
      cellscatra_subproblem = Teuchos::rcp(new SSI::SSI_Part2WC_PROTRUSIONFORMATION(comm, problem->CellMigrationParams(), problem->CellMigrationParams().sublist("SCALAR TRANSPORT"), problem->CellMigrationParams().sublist("STRUCTURAL DYNAMIC"), "cell", "cellscatra"));
      break;
    default:
      dserror("unknown coupling algorithm for SSI! Only ssi_IterStagg valid. Fix your *.dat file.");
      break;
    }

    // set pointer to structure inside SSI subproblem
    cellstructure = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapperImmersed>(cellscatra_subproblem->StructureField());

    if(cellstructure==Teuchos::null)
      dserror("dynamic cast from Structure to FSIStructureWrapperImmersed failed");

    if(comm.MyPID()==0)
      std::cout<<"\nCreated Field Cell Structure with intracellular signaling capabilitiy...\n \n"<<std::endl;
  }
  else if(not ssi_cell) // cell has no intracellular biochemistry -> just create usual fsi immersed structure
  {
    cellstructure = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapperImmersed>(Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(problem->CellMigrationParams(),
                                                                                                                                     (problem->CellMigrationParams()).sublist("STRUCTURAL DYNAMIC"),
                                                                                                                                     problem->GetDis("cell")))->StructureField() );

    if(cellstructure==Teuchos::null)
      dserror("dynamic cast from Structure to FSIStructureWrapperImmersed failed");

    if(comm.MyPID()==0)
      std::cout<<"\n Created Field Cell Structure without intracellular signaling capabilitiy... \n \n"<<std::endl;
  }
  else
  {
      std::cout<<"Skipped Construction of CellScatra Subproblem.\n"
              "The current parameter combination is not supported.\n"
              "Check your .dat file !"<<std::endl;;
  }

  // create instance of poroelast subproblem
  Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem = POROELAST::UTILS::CreatePoroScatraAlgorithm(problem->CellMigrationParams(),comm);

  // setup of poro monolithic algorithm
  poroscatra_subproblem->SetupSystem();
  if(comm.MyPID()==0)
    std::cout<<" Created Field PoroScatra ... \n"<<std::endl;



  //////////////////////////////////////////////
  // setup the immersed discretization
  //////////////////////////////////////////////
  std::map<int,LINALG::Matrix<3,1> > currpositions_cell; //!< map of vectors for search tree containing current structural positions

  // construct 3D search tree for cell domain
  Teuchos::RCP<GEO::SearchTree> cell_SearchTree = Teuchos::rcp(new GEO::SearchTree(5));

  // ghost structure on each proc (for search algorithm)
  if(comm.NumProc() > 1 and simtype!=INPAR::CELL::sim_type_pureContraction)
  {
    // fill complete inside
    CreateGhosting(problem->GetDis("cell"));
  }
  else
  {
    // fill complete to incorporate changes due to ghosting and build geometries
    problem->GetDis("cell")->FillComplete();
  }

  // find positions of the immersed discretization on proc
  std::map<int,LINALG::Matrix<3,1> > my_currpositions_cell;
  for (int lid = 0; lid < problem->GetDis("cell")->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = problem->GetDis("cell")->lRowNode(lid);
    LINALG::Matrix<3,1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    my_currpositions_cell[node->Id()] = currpos;
  }
  // Communicate local currpositions:
  // map with current structural positions should be same on all procs
  // to make use of the advantages of ghosting the structure redundantly
  // on all procs.
  int procs[comm.NumProc()];
  for(int i=0;i<comm.NumProc();i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_cell,currpositions_cell,comm.NumProc(),&procs[0],comm);

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox2 = GEO::getXAABBofDis(*(problem->GetDis("cell")),currpositions_cell);
  cell_SearchTree->initializeTree(rootBox2,*(problem->GetDis("cell")),GEO::TreeType(GEO::OCTTREE));

  if(comm.MyPID()==0)
    std::cout<<"\n Build Cell SearchTree ... "<<std::endl;



  //////////////////////////////////////////////
  // setup background discretization
  //////////////////////////////////////////////
  std::map<int,LINALG::Matrix<3,1> > currpositions_ECM;  //!< pointer to map of vectors for search tree containing current fluid positions

  // construct 3D search tree for fluid domain
  Teuchos::RCP<GEO::SearchTree> fluid_SearchTree = Teuchos::rcp(new GEO::SearchTree(5));

  // find positions of the background fluid discretization
  for (int lid = 0; lid < problem->GetDis("porofluid")->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = problem->GetDis("porofluid")->lColNode(lid);
    LINALG::Matrix<3,1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    currpositions_ECM[node->Id()] = currpos;
  }

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*(problem->GetDis("porofluid")),currpositions_ECM);
  fluid_SearchTree->initializeTree(rootBox,*(problem->GetDis("porofluid")),GEO::TreeType(GEO::OCTTREE));

  if(comm.MyPID()==0)
    std::cout<<"\n Build Fluid/ECM SearchTree ... "<<std::endl;



  // write pointers to previously created objects into ParameterList in order to make them
  // accessible in subproblems
  params.set<Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> >("RCPToCellStructure",cellstructure);
  params.set<Teuchos::RCP<POROELAST::PoroScatraBase> >("RCPToPoroScatra",poroscatra_subproblem);
  params.set<Teuchos::RCP<GEO::SearchTree> >("RCPToCellSearchTree",cell_SearchTree);
  params.set<Teuchos::RCP<GEO::SearchTree> >("RCPToFluidSearchTree",fluid_SearchTree);
  params.set<std::map<int,LINALG::Matrix<3,1> >* >("PointerToCurrentPositionsCell",&currpositions_cell);
  params.set<std::map<int,LINALG::Matrix<3,1> >* >("PointerToCurrentPositionsECM",&currpositions_ECM);
  params.set<Teuchos::RCP<SSI::SSI_Part2WC_PROTRUSIONFORMATION> >("RCPToCellScatra",cellscatra_subproblem);

  //////////////////////////////////////////////
  // query simulation type
  //////////////////////////////////////////////
  if(simtype==INPAR::CELL::sim_type_pureFSI)
  {
    Teuchos::RCP<IMMERSED::ImmersedPartitionedFlowCellInteraction> algo =
        Teuchos::rcp(new IMMERSED::ImmersedPartitionedFlowCellInteraction(params,comm));

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      algo->ReadRestart(restart);
    }

    algo->Timeloop(algo);

    if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
    {
      Teuchos::TimeMonitor::summarize();
      Teuchos::TimeMonitor::zeroOutTimers();
    }

    // create result tests for single fields
    DRT::Problem::Instance()->AddFieldTest(algo->CellField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(algo->PoroField()->FluidField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(algo->PoroField()->StructureField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(algo->PoroField()->ScaTraFieldBase()->CreateScaTraFieldTest());

    // do the actual testing
    DRT::Problem::Instance()->TestAll(comm);

  }// sim_type_pureFSI
  else if (simtype==INPAR::CELL::sim_type_pureAdhesion)
  {
    params.set<bool>("IsPureAdhesionSimulation", true);

    Teuchos::RCP<IMMERSED::ImmersedPartitionedAdhesionTraction> algo =
        Teuchos::rcp(new IMMERSED::ImmersedPartitionedAdhesionTraction(params,comm));

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      algo->ReadRestart(restart);
    }

    algo->Timeloop(algo);

    if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
    {
      Teuchos::TimeMonitor::summarize();
      Teuchos::TimeMonitor::zeroOutTimers();
    }

    // create result tests for single fields
    DRT::Problem::Instance()->AddFieldTest(cellstructure->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->FluidField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->StructureField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->ScaTraFieldBase()->CreateScaTraFieldTest());

    // do the actual testing
    DRT::Problem::Instance()->TestAll(comm);

  }// sim_type_pureAdhesion
  else if (simtype==INPAR::CELL::sim_type_pureConfinement)
  {
    params.set<bool>("IsPureConfinementSimulation", true);

    Teuchos::RCP<IMMERSED::ImmersedPartitionedConfineCell> algo =
        Teuchos::rcp(new IMMERSED::ImmersedPartitionedConfineCell(params,comm));

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      algo->ReadRestart(restart);
    }

    algo->Timeloop(algo);

    if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
    {
      Teuchos::TimeMonitor::summarize();
      Teuchos::TimeMonitor::zeroOutTimers();
    }

    // create result tests for single fields
    DRT::Problem::Instance()->AddFieldTest(cellstructure->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->FluidField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->StructureField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->ScaTraFieldBase()->CreateScaTraFieldTest());

    // do the actual testing
    DRT::Problem::Instance()->TestAll(comm);

  }// sim_type_pureCompression
  else if (simtype==INPAR::CELL::sim_type_pureProtrusionFormation)
  {
    params.set<bool>("IsPureProtrusionFormation", true);

    Teuchos::RCP<IMMERSED::ImmersedPartitionedProtrusionFormation> algo =
        Teuchos::rcp(new IMMERSED::ImmersedPartitionedProtrusionFormation(params,comm));

    // ghost cellscatra on each proc (for search algorithm)
    if(comm.NumProc() > 1)
    {
      // fill complete inside
      CreateGhosting(problem->GetDis("cellscatra"));
    }

    // handle restart
    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      algo->ReadRestart(restart);
    }

    // run the problem
    algo->Timeloop(algo);

    if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
    {
      Teuchos::TimeMonitor::summarize();
      Teuchos::TimeMonitor::zeroOutTimers();
    }

    // create result tests for single fields
    DRT::Problem::Instance()->AddFieldTest(cellstructure->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->FluidField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->StructureField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(poroscatra_subproblem->ScaTraFieldBase()->CreateScaTraFieldTest());

    // do the actual testing
    DRT::Problem::Instance()->TestAll(comm);

  }// sim_type_pureProtrusionFormation
  else if (simtype==INPAR::CELL::sim_type_pureContraction)
  {
    if(comm.MyPID()==0)
      std::cout<<"\n Simulation Type : Pure Biochemo-Mechano Coupled Cell Contraction  \n"<<std::endl;

    params.set<bool>("IsPureContraction", true);

    Teuchos::RCP<SSI::SSI_Part2WC_BIOCHEMOMECHANO> algo =
        Teuchos::rcp(new SSI::SSI_Part2WC_BIOCHEMOMECHANO(comm, params, problem->CellMigrationParams(), problem->ScalarTransportDynamicParams(), problem->StructuralDynamicParams(), "cell", "cellscatra"));

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      algo->ReadRestart(restart);
    }

    // run the problem
    algo->Timeloop();

    if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
    {
      Teuchos::TimeMonitor::summarize();
      Teuchos::TimeMonitor::zeroOutTimers();
    }

    // perform the result test
    algo->TestResults(comm);

  }// sim_type_pureContraction
  else if (simtype==INPAR::CELL::sim_type_pureEndoExocytosis)
  {
    scatra_dyn(0);
  }
  else if (simtype==INPAR::CELL::sim_type_multiphysics)
  {
    // here the global algo of the fully coupled cell migration will be implemented
  }
  else
    dserror("Unknown SIMTYPE set in ---CELL DYNAMIC section (default=pureFSI). Fix your .dat file.");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CreateGhosting(const Teuchos::RCP<DRT::Discretization> distobeghosted)
{
  if(distobeghosted->Comm().MyPID() == 0)
  {
    std::cout<<"################################################################################################"<<std::endl;
    std::cout<<"###   Ghost discretization "<<distobeghosted->Name()<<" redundantly on all procs ... "<<std::endl;
    std::cout<<"################################################################################################"<<std::endl;
  }

  std::vector<int> allproc(distobeghosted->Comm().NumProc());
  for (int i=0; i<distobeghosted->Comm().NumProc(); ++i) allproc[i] = i;

  // fill my own row node ids
  const Epetra_Map* noderowmap = distobeghosted->NodeRowMap();
  std::vector<int> sdata;
  for (int i=0; i<noderowmap->NumMyElements(); ++i)
  {
    int gid = noderowmap->GID(i);
    DRT::Node* node = distobeghosted->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    sdata.push_back(gid);
  }

  // gather all master row node gids redundantly
  std::vector<int> rdata;
  LINALG::Gather<int>(sdata,rdata,(int)allproc.size(),&allproc[0],distobeghosted->Comm());

  // build new node column map (on ALL processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,distobeghosted->Comm()));
  sdata.clear();
  rdata.clear();

  // fill my own row element ids
  const Epetra_Map* elerowmap  = distobeghosted->ElementRowMap();
  sdata.resize(0);
  for (int i=0; i<elerowmap->NumMyElements(); ++i)
  {
    int gid = elerowmap->GID(i);
    DRT::Element* ele = distobeghosted->gElement(gid);
    if (!ele) dserror("ERROR: Cannot find element with gid %",gid);
    sdata.push_back(gid);
  }

  // gather all gids of elements redundantly
  rdata.resize(0);
  LINALG::Gather<int>(sdata,rdata,(int)allproc.size(),&allproc[0],distobeghosted->Comm());

  // build new element column map (on ALL processors)
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,distobeghosted->Comm()));
  sdata.clear();
  rdata.clear();
  allproc.clear();

  // redistribute the discretization of the interface according to the
  // new node / element column layout (i.e. master = full overlap)
  distobeghosted->ExportColumnNodes(*newnodecolmap);
  distobeghosted->ExportColumnElements(*newelecolmap);

  // finalize the changes made to the discretisation
  distobeghosted->FillComplete();

#ifdef DEBUG
  int nummycolnodes = newnodecolmap->NumMyElements();
  int nummycolelements = newelecolmap->NumMyElements();
  int sizelist[distobeghosted->Comm().NumProc()];
  distobeghosted->Comm().GatherAll(&nummycolnodes,&sizelist[0],1);
  std::cout<<"PROC "<<distobeghosted->Comm().MyPID()<<" : "<<nummycolnodes<<" colnodes"<<std::endl;
  distobeghosted->Comm().Barrier();
  std::cout<<"PROC "<<distobeghosted->Comm().MyPID()<<" : "<<nummycolelements<<" colelements"<<std::endl;
  distobeghosted->Comm().Barrier();
  std::cout<<"PROC "<<distobeghosted->Comm().MyPID()<<" : "<<distobeghosted->lColElement(0)->Nodes()[0]->Id()<<" first ID of first node of first colele"<<std::endl;
  distobeghosted->Comm().Barrier(); // wait for procs
  for(int k=1;k<distobeghosted->Comm().NumProc();++k)
  {
    if(sizelist[k-1]!=nummycolnodes)
      dserror("Since whole dis is ghosted every processor should have the same number of colnodes. This is not the case! Fix this!");
  }
#endif

}

