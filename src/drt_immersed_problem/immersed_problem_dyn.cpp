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
#include "immersed_partitioned_cellcontraction.H"
#include "immersed_partitioned_adhesion_traction.H"
#include "ssi_partitioned_2wc_protrusionformation.H"
#include "immersed_partitioned_fsi_dirichletneumann.H"
#include "immersed_partitioned_protrusion_formation.H"
#include "immersed_partitioned_flow_cell_interaction.H"
#include "immersed_partitioned_fsi_dirichletneumann_ale.H"
#include "str_model_evaluator_multiphysics_cellmigration.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_inpar/inpar_immersed.H"
#include "../drt_inpar/inpar_cell.H"
#include "../drt_inpar/inpar_ssi.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_multiphysicswrapper_cellmigration.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_ale/ale_utils_clonestrategy.H"
#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_poroelast/poroelast_utils_setup.H"
#include "../drt_poroelast/poroelast_monolithic.H"
#include "../drt_poroelast/poro_scatra_part_2wc.H"
#include "../drt_scatra/scatra_dyn.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_heterogeneous_reaction_strategy.H"



void immersed_problem_drt()
{
  ///////////////////////////////////////////////////////////////////////
  // General declarations and variables
  ///////////////////////////////////////////////////////////////////////
  // declare general ParameterList that can be handed into Algorithm
  Teuchos::ParameterList params;
  // get pointer to global problem
  DRT::Problem* problem = DRT::Problem::Instance();
  // get communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  ///////////////////////////////////////////////////////////////////////
  // Get Parameters
  ///////////////////////////////////////////////////////////////////////
  // get parameterlist for immersed method
  const Teuchos::ParameterList& immersedmethodparams = problem->ImmersedMethodParams();
  // choose algorithm
  int coupling = DRT::INPUT::IntegralValue<int>(immersedmethodparams,"COUPALGO");
  int scheme   = DRT::INPUT::IntegralValue<int>(immersedmethodparams,"SCHEME");


  ///////////////////////////////////////////////////////////////////////
  // Query algorithms
  ///////////////////////////////////////////////////////////////////////
  switch (coupling)
  {
  case INPAR::IMMERSED::partitioned:
  {
    switch(DRT::Problem::Instance()->ProblemType())
    {
    case prb_immersed_fsi:
    {
      // fill discretizations
      problem->GetDis("structure")->FillComplete(false,false,false);
      problem->GetDis("fluid")    ->FillComplete(false,false,false);

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

      Teuchos::RCP<IMMERSED::ImmersedPartitionedFSIDirichletNeumann> algo = Teuchos::null;
      if(scheme == INPAR::IMMERSED::dirichletneumann)
        algo = Teuchos::rcp(new IMMERSED::ImmersedPartitionedFSIDirichletNeumann(comm));
      else
      {
        algo = Teuchos::null;
        dserror("unknown coupling scheme");
      }

      // set fluid action type needed in function CalcArtificialVelocity()
      params.set<std::string>("fluid_action","interpolate_velocity_to_given_point_immersed");

      // init algo
      algo->Init(params);

      // ghost structure redundantly on all procs
      DRT::UTILS::GhostDiscretizationOnAllProcs(problem->GetDis("structure"));

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
      // fill discretizations
      problem->GetDis("structure")->FillComplete(false,false,false);
      problem->GetDis("fluid")    ->FillComplete(false,false,false);
      problem->GetDis("ale")      ->FillComplete(false,false,false);

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

      // get discretizations
      Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
      Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

      // create ale elements if the ale discretization is empty
      if (aledis->NumGlobalNodes()==0)
      {
        DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
        aledis->FillComplete();
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

      // set fluid action type needed in function CalcArtificialVelocity()
      params.set<std::string>("fluid_action","interpolate_velocity_to_given_point_immersed");

      // init algo
      algo->Init(params);

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
  ///////////////////////////////////////////////////////////////////////
  // General declarations and variables
  ///////////////////////////////////////////////////////////////////////
  // declare general ParameterList that can be handed into Algorithm
  Teuchos::ParameterList params;
  // declare pointer to field cell structure
  Teuchos::RCP< ::ADAPTER::MultiphysicsStructureWrapperCellMigration>
      cellstructure = Teuchos::null;
  // declare model evaluator
  Teuchos::RCP<STR::MODELEVALUATOR::CellMigration>
      multiphysics_model_evaluator_cellmigration = Teuchos::null;
  // declare pointer to cell subproblem (structure-scatra interaction)
  Teuchos::RCP<SSI::SSI_Part2WC> cellscatra_subproblem = Teuchos::null;

  // get pointer to global problem
  DRT::Problem* problem = DRT::Problem::Instance();
  // get communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  ///////////////////////////////////////////////////////////////////////
  // Get Parameters
  ///////////////////////////////////////////////////////////////////////
  // get parameterlist for immersed method
  const Teuchos::ParameterList& immersedmethodparams = problem->ImmersedMethodParams();
  // get parameter list for cell migration
  const Teuchos::ParameterList& cellmigrationparams = problem->CellMigrationParams();
  // get parameter list for cell structure
  const Teuchos::ParameterList& cellstructureparams = cellmigrationparams.sublist("STRUCTURAL DYNAMIC");
  // extract the simulation type
  int simtype = DRT::INPUT::IntegralValue<int>(cellmigrationparams,"SIMTYPE");
  // decides if cell is constructed as pure structure or SSI algorithm
  bool ssi_cell = DRT::INPUT::IntegralValue<int>(problem->CellMigrationParams(),"SSI_CELL");

  // use PoroelastImmersedCloneStrategy to build FluidPoroImmersed elements
  POROELAST::UTILS::SetupPoroScatraDiscretizations
  <POROELAST::UTILS::PoroelastImmersedCloneStrategy,POROELAST::UTILS::PoroScatraCloneStrategy>();

  ///////////////////////////////////////////////////////////////////////
  // Safety check
  ///////////////////////////////////////////////////////////////////////
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

  ///////////////////////////////////////////////////////////////////////
  // Initial FillComplete() just creates basic connectivity
  ///////////////////////////////////////////////////////////////////////
  problem->GetDis("cell")->FillComplete(true,true,true);
  problem->GetDis("cellscatra")->FillComplete(true,true,true);



  ///////////////////////////////////////////////////////////////////////
  // Construct the cell structure
  ///////////////////////////////////////////////////////////////////////
  // construct base algorithm
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew>
      struct_adapterbase_ptr = ADAPTER::STR::BuildStructureAlgorithm(cellstructureparams);
  // initialize structure base algorithm
  struct_adapterbase_ptr->Init(
      cellmigrationparams,
      const_cast<Teuchos::ParameterList&>(cellstructureparams),
      problem->GetDis("cell"));



  ///////////////////////////////////////////////////////////////////////
  // Construct model evaluators for cell migration
  ///////////////////////////////////////////////////////////////////////
  // construct model evaluator for cell migration
  multiphysics_model_evaluator_cellmigration =
      Teuchos::rcp(new STR::MODELEVALUATOR::CellMigration());

  // register model evaluator in adapter structure base object
  struct_adapterbase_ptr -> RegisterModelEvaluator(
      "Partitioned Coupling Model",
      multiphysics_model_evaluator_cellmigration );



  ///////////////////////////////////////////////////////////////////////
  // Construct cell
  ///////////////////////////////////////////////////////////////////////
  // check if cell is supposed to be constructed as SSI algo,
  // i.e., with intracellular biochchemical signaling
  if(ssi_cell)
  {
    // we need ssi with ale formulation, since our cell is a deforming domain
    bool isale = true;

    // get coupling algorithm from ---SSI CONTROL section
    const INPAR::SSI::SolutionSchemeOverFields coupling
    = DRT::INPUT::IntegralValue<INPAR::SSI::SolutionSchemeOverFields>(problem->SSIControlParams(),"COUPALGO");

    // construct specific SSI type
    switch(coupling)
    {
    case INPAR::SSI::ssi_IterStagg:
      if(simtype==INPAR::CELL::sim_type_pureProtrusionFormation)
        cellscatra_subproblem = Teuchos::rcp(new SSI::SSI_Part2WC_PROTRUSIONFORMATION(comm,problem->CellMigrationParams()));
      else if (simtype==INPAR::CELL::sim_type_pureAdhesion)
        cellscatra_subproblem = Teuchos::rcp(new SSI::SSI_Part2WC(comm,problem->CellMigrationParams()));
      else if (simtype==INPAR::CELL::sim_type_pureContraction)
      {
        cellscatra_subproblem = Teuchos::rcp(new SSI::SSI_Part2WC_BIOCHEMOMECHANO(comm,problem->CellMigrationParams()));
        params.set<bool>("IsPureContraction", true);
      }
        break;
    default:
      dserror("unknown coupling algorithm for SSI! Only ssi_IterStagg valid. Fix your *.dat file.");
      break;
    }

    // set the pointer to the structure base algo in SSI object
    cellscatra_subproblem->SetStructureAdapterBase(struct_adapterbase_ptr);

    // It is time to call Init().
    // "ale" dis is cloned and filled inside.
    // SSI coupling object is constructed inside.
    // Returns true if redistribution is required.
      int redistribute =
      cellscatra_subproblem->Init(
          comm,
          problem->CellMigrationParams(),
          problem->CellMigrationParams().sublist("SCALAR TRANSPORT"),
          problem->CellMigrationParams().sublist("STRUCTURAL DYNAMIC"),
          "cell",
          "cellscatra",
          isale);

    // Redistribute discretizations if necessary
    if(redistribute)
    {
      if( redistribute== (int)SSI::match)
      {
        problem->GetDis("cell")->FillComplete(true,true,true);
        problem->GetDis("cellscatra")->FillComplete(true,true,true);

        if(simtype==INPAR::CELL::sim_type_pureProtrusionFormation)
          problem->GetDis("ale")->FillComplete(true,true,true);

        // first we bin the scatra discretization
        std::vector<Teuchos::RCP<DRT::Discretization> > dis;
        dis.push_back(problem->GetDis("cellscatra"));
        DRT::UTILS::RedistributeDiscretizationsByBinning(dis,false);

        DRT::UTILS::MatchElementDistributionOfMatchingConditionedElements(
            *problem->GetDis("cellscatra"),
            *problem->GetDis("cellscatra"),
            "ScatraHeteroReactionMaster",
            "ScatraHeteroReactionSlave" );

        // now we redistribute the structure dis to match the scatra dis
        DRT::UTILS::MatchElementDistributionOfMatchingDiscretizations(
            *problem->GetDis("cellscatra"),
            *problem->GetDis("cell") );

        // now we redistribute the ale dis to match the scatra dis
        if(simtype==INPAR::CELL::sim_type_pureProtrusionFormation)
          DRT::UTILS::MatchElementDistributionOfMatchingDiscretizations(
              *problem->GetDis("cell"),
              *problem->GetDis("ale") );
      }
      else
        dserror("Only matching discretizatiosn are supported.");

    } // end redistribution


    // ghost cell discretizations on each proc (for search algorithm)
    if(comm.NumProc() > 1)
    {
      // fill complete inside
      DRT::UTILS::GhostDiscretizationOnAllProcs(problem->GetDis("cell"));
      DRT::UTILS::GhostDiscretizationOnAllProcs(problem->GetDis("cellscatra"));
      if(simtype==INPAR::CELL::sim_type_pureProtrusionFormation)
        DRT::UTILS::GhostDiscretizationOnAllProcs(problem->GetDis("ale"));
    } // end ghosting


    // now we call the final fill complete on our discretizations.
    // FillComplete for ale dis is called deeper in the code.
    problem->GetDis("cell")->FillComplete(true,true,true);
    problem->GetDis("cellscatra")->FillComplete(true,true,true);
    problem->GetDis("ale")->FillComplete(true,true,true);



    // parallel redistribution is finished. Let us call Setup()
    // here all state vectors are constructed.
    // now we are sure, that all maps fit the actual distribution.
    //(wrapper is created inside)
    struct_adapterbase_ptr -> Setup();

    // extract the problem specific wrapper
    cellstructure =
        Teuchos::rcp_dynamic_cast<ADAPTER::MultiphysicsStructureWrapperCellMigration>(
            struct_adapterbase_ptr->StructureField(),true);

    // set pointer to model evaluator in ssi specific wrapper
    cellstructure->GetSSIStructureWrapperPtr()->
        SetModelEvaluatorPtr(
            Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::PartitionedSSI>(
            multiphysics_model_evaluator_cellmigration->GetModelEvaluatorFromMap(STR::MODELEVALUATOR::mt_ssi),
            true) );

    // set the ssi specific wrapper in SSI object
    cellscatra_subproblem -> SetStructureWrapper(cellstructure->GetSSIStructureWrapperPtr());

    // set the struct-ale specific wrapper in SSI object for protruison formation
    if(simtype==INPAR::CELL::sim_type_pureProtrusionFormation)
      Teuchos::rcp_dynamic_cast<SSI::SSI_Part2WC_PROTRUSIONFORMATION>(cellscatra_subproblem)
          -> SetStructAleWrapper(cellstructure->GetStructAleStructureWrapperPtr());

    // now setup cellscatra (ssi) subproblem
    cellscatra_subproblem->Setup();

    if(comm.MyPID()==0)
      std::cout<<"\nCreated Field Cell Structure with intracellular signaling capabilitiy...\n \n"<<std::endl;

  }
  else // cell has no intracellular biochemistry -> just create usual structure
  {
    // ghosting
    if(comm.NumProc() > 1)
    {
      DRT::UTILS::GhostDiscretizationOnAllProcs(problem->GetDis("cell"));
      problem->GetDis("cell")->FillComplete(true,true,true);
    }

    // setup structure (wrapper is created)
    struct_adapterbase_ptr -> Setup();

    // extract the problem specific wrapper
    cellstructure =
        Teuchos::rcp_dynamic_cast<ADAPTER::MultiphysicsStructureWrapperCellMigration>(
            struct_adapterbase_ptr->StructureField(),true);

    if(comm.MyPID()==0)
      std::cout<<"\n Created Field Cell Structure without intracellular signaling capabilitiy... \n \n"<<std::endl;

  } // finish construction of cell

  ///////////////////////////////////////////////////////////////////////
  // Safety check
  ///////////////////////////////////////////////////////////////////////
  if(cellstructure==Teuchos::null)
    dserror("dynamic cast from Structure to MultiphysicsStructureWrapperCellMigration failed");

  // set pointer to model evaluator in ssi specific wrapper
  cellstructure->GetFSIStructureWrapperPtr()->
      SetModelEvaluatorPtr(
          Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::PartitionedFSI>(
          multiphysics_model_evaluator_cellmigration->GetModelEvaluatorFromMap(STR::MODELEVALUATOR::mt_fsi),
          true) );



  ///////////////////////////////////////////////////////////////////////
  // Construct the poro-scatra subproblem for the extracellular matrix
  ///////////////////////////////////////////////////////////////////////
  // create instance of poroelast subproblem
  Teuchos::RCP<POROELAST::PoroScatraBase>
      poroscatra_subproblem = POROELAST::UTILS::CreatePoroScatraAlgorithm(
          problem->CellMigrationParams(),
          comm);

  // setup of poro monolithic algorithm
  poroscatra_subproblem->SetupSystem();
  if(comm.MyPID()==0)
    std::cout<<" Created Field PoroScatra ... \n"<<std::endl;



  ///////////////////////////////////////////////////////////////////////
  // Setup the immersed discretization
  ///////////////////////////////////////////////////////////////////////
  // map of vectors for search tree containing current structural positions
  std::map<int,LINALG::Matrix<3,1> > currpositions_cell;

  // construct 3D search tree for cell domain
  Teuchos::RCP<GEO::SearchTree> cell_SearchTree = Teuchos::rcp(new GEO::SearchTree(5));

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



  ///////////////////////////////////////////////////////////////////////
  // Setup background discretization
  ///////////////////////////////////////////////////////////////////////
  // pointer to map of vectors for search tree containing current ECM positions
  std::map<int,LINALG::Matrix<3,1> > currpositions_ECM;


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
    std::cout<<"\n Build Fluid/ECM SearchTree ... \n"<<std::endl;



  ///////////////////////////////////////////////////////////////////////
  // provide pointers to previously created objects to ParameterList
  // in order to make them accessible in subproblems.
  ///////////////////////////////////////////////////////////////////////
  params.set<Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration> >("RCPToCellStructure",cellstructure);
  params.set<Teuchos::RCP<POROELAST::PoroScatraBase> >("RCPToPoroScatra",poroscatra_subproblem);
  params.set<Teuchos::RCP<GEO::SearchTree> >("RCPToCellSearchTree",cell_SearchTree);
  params.set<Teuchos::RCP<GEO::SearchTree> >("RCPToFluidSearchTree",fluid_SearchTree);
  params.set<std::map<int,LINALG::Matrix<3,1> >* >("PointerToCurrentPositionsCell",&currpositions_cell);
  params.set<std::map<int,LINALG::Matrix<3,1> >* >("PointerToCurrentPositionsECM",&currpositions_ECM);
  params.set<Teuchos::RCP<SSI::SSI_Part2WC> >("RCPToCellScatra",cellscatra_subproblem);




  ///////////////////////////////////////////////////////////////////////
  // Query simulation type
  ///////////////////////////////////////////////////////////////////////
  if(simtype==INPAR::CELL::sim_type_pureFSI)
  {
    Teuchos::RCP<IMMERSED::ImmersedPartitionedFlowCellInteraction> algo =
        Teuchos::rcp(new IMMERSED::ImmersedPartitionedFlowCellInteraction(params,comm));

    // init algo
    algo->Init(params);

    // setup algo
    algo->Setup();

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

    // init algo
    algo->Init(params);

    // setup algo
    algo->Setup();

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

    // init algo
    algo->Init(params);

    // setup algo
    algo->Setup();

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

    // init algo
    algo->Init(params);

    // setup algo
    algo->Setup();

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
    params.set<bool>("IsPureContractionSimulation", true);

    Teuchos::RCP<IMMERSED::ImmersedPartitionedCellContraction> algo =
        Teuchos::rcp(new IMMERSED::ImmersedPartitionedCellContraction(params,comm));

    // init algo
    algo->Init(params);

    // setup algo
    algo->Setup();

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      algo->ReadRestart(restart);
    }

    // run the problem
    algo -> Timeloop(algo);

    if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
    {
      Teuchos::TimeMonitor::summarize();
      Teuchos::TimeMonitor::zeroOutTimers();
    }

    // create result tests for single fields
    DRT::Problem::Instance()->AddFieldTest(cellscatra_subproblem->ScaTraField()->CreateScaTraFieldTest());

    // do the actual testing
    DRT::Problem::Instance()->TestAll(comm);

  }// sim_type_pureContraction

  else if (simtype==INPAR::CELL::sim_type_pureEndoExocytosis)
  {
    scatra_dyn(0);
  }// sim_type_pureEndoExocytosis

  else if (simtype==INPAR::CELL::sim_type_multiphysics)
  {
    // here the global algo of the fully coupled cell migration will be implemented
  }// sim_type_multiphysics

  else
    dserror("Unknown SIMTYPE set in ---CELL DYNAMIC section (default=pureFSI). Fix your .dat file.");

  return;
}

