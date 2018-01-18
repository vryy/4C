/*!------------------------------------------------------------------------------------------------*
 \file ssi_base.cpp

 \brief base class for all scalar structure algorithms

 \level 1

   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de

 *------------------------------------------------------------------------------------------------*/

#include "ssi_base.H"

#include "ssi_partitioned.H"
#include "ssi_coupling.H"
#include "ssi_resulttest.H"
#include "ssi_utils.H"
#include "ssi_str_model_evaluator_partitioned.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
//for cloning
#include "../drt_lib/drt_utils_createdis.H"
#include "ssi_clonestrategy.H"

#include "../drt_inpar/inpar_s2i.H"
#include "../drt_inpar/inpar_volmortar.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra_ele/scatra_ele.H"

#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Base::SSI_Base(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams
    ) :
    AlgorithmBase(comm, globaltimeparams),
    icoup_structure_(Teuchos::null),
    iter_(0),
    map_structure_condensed_(Teuchos::null),
    maps_structure_(Teuchos::null),
    struct_adapterbase_ptr_(Teuchos::null),
    structure_(Teuchos::null),
    scatra_(Teuchos::null),
    ssicoupling_(Teuchos::null),
    use_old_structure_(false), // todo temporary flag
    zeros_(Teuchos::null),
    fieldcoupling_(DRT::INPUT::IntegralValue<INPAR::SSI::FieldCoupling>(DRT::Problem::Instance()->SSIControlParams(),"FIELDCOUPLING")),
    issetup_(false),
    isinit_(false)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call the setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Special Init in case adapter needs to be set from outside            |
 |                                                          rauch 11/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::SetStructureAdapterBase(
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> struct_adapterbase_ptr)
{
  struct_adapterbase_ptr_ = struct_adapterbase_ptr;
  return;
}

/*----------------------------------------------------------------------*
 | Special init in case adapter needs to be set from outside            |
 |                                                          rauch 11/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::SetStructureWrapper(
    Teuchos::RCP<ADAPTER::Structure> struct_wrapper)
{
  structure_ =
      Teuchos::rcp_dynamic_cast<ADAPTER::SSIStructureWrapper>(struct_wrapper);
  return;
}

/*----------------------------------------------------------------------*
 | Init this class                                          rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Base::Init(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname,
    bool isAle
    )
{
  // reset the setup flag
  SetIsSetup(false);

  // get the global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatraparams.get<int>("LINEAR_SOLVER");

  //3.- Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis(struct_disname);

  // !!! TIME PARAMETER HANDLING !!!
  // Determine which time params to use to build the single fields.
  // In case of time stepping, time params have to be read from single field sections.
  // In case of equal timestep size for all fields the time params are controlled solely
  // by the problem section (e.g. ---SSI DYNAMIC or ---CELL DYNAMIC).
  //
  // access the ssi dynamic params
  const Teuchos::ParameterList* structtimeparams = &globaltimeparams;
  const Teuchos::ParameterList* scatratimeparams = &globaltimeparams;
  if(DRT::INPUT::IntegralValue<int>(globaltimeparams,"DIFFTIMESTEPSIZE"))
  {
    structtimeparams = &structparams;
    scatratimeparams = &scatraparams;
  }

  // we do not construct a structure here, in case it
  // was built externally and handed into this object.
  if (struct_adapterbase_ptr_ == Teuchos::null)
  {
    // create structural field
    // todo FIX THIS !!!!
    // Decide whether to use old structural time integration or new structural time integration.
    // This should be removed as soon as possible!
    // Also all structural elements need to be adapted first!
    // Then, we can switch the remaining ssi tests using the old time integration to the new one,
    // i.e.: todo
    //
    // build structure
    if(structparams.get<std::string>("INT_STRATEGY")=="Standard")
    {
      struct_adapterbase_ptr_= ADAPTER::STR::BuildStructureAlgorithm(structparams);

      // initialize structure base algorithm
      struct_adapterbase_ptr_->Init(
          *structtimeparams,
          const_cast<Teuchos::ParameterList&>(structparams),
          structdis);
    }
    else if (structparams.get<std::string>("INT_STRATEGY")=="Old") // todo this is the part that should be removed !
    {
      if (comm.MyPID()==0)
        std::cout<<"\n"<<
        " USING OLD STRUCTURAL TIME INTEGRATION FOR STRUCTURE-SCATRA-INTERACTION!\n"
        " FIX THIS! THIS IS ONLY SUPPOSED TO BE TEMPORARY!"
        "\n"<<std::endl;

      use_old_structure_ = true;

      Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
          Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(*structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis));
      structure_ = Teuchos::rcp_dynamic_cast< ::ADAPTER::SSIStructureWrapper>(structure->StructureField(),true);
      if(structure_ == Teuchos::null)
        dserror("cast from ADAPTER::Structure to ADAPTER::SSIStructureWrapper failed");
    }
    else
      dserror("Unknown time integration requested!\n"
          "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!\n"
          "If you want to use yet unsupported elements or algorithms,\n"
          "set INT_STRATEGY to Old in ---STRUCUTRAL DYNAMIC section!");

  } // if structure_ not set from outside

  // create scatra field
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

  // do discretization specific setup (e.g. clone discr. scatra from structure)
  InitDiscretizations(comm,struct_disname,scatra_disname);

  // initialize scatra base algorithm.
  // scatra time integrator constructed and initialized inside.
  // mesh is written inside. cloning must happen before!
  scatra_->Init(
      *scatratimeparams,
      scatraparams,
      problem->SolverParams(linsolvernumber),
      scatra_disname,
      isAle);

  int redistribute = InitFieldCoupling(comm,struct_disname,scatra_disname);

  // set isinit_ flag true
  SetIsInit(true);

  return redistribute;
}

/*----------------------------------------------------------------------*
 | Setup this class                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::Setup()
{
  // check initialization
  CheckIsInit();

  // set up helper class for field coupling
  ssicoupling_->Setup();

  // set up scalar transport field
  scatra_->ScaTraField()->Setup();

  // only relevant for new structural time integration
  // only if adapter base has not already been set up outside
  if(not use_old_structure_ and not struct_adapterbase_ptr_->IsSetup())
  {
    // set up structural model evaluator
    SetupModelEvaluator();

    // pass initial scalar field to structural discretization to correctly compute initial accelerations
    if(DRT::INPUT::IntegralValue<INPAR::SSI::SolutionSchemeOverFields>(DRT::Problem::Instance()->SSIControlParams(),"COUPALGO") != INPAR::SSI::ssi_OneWay_SolidToScatra)
      ssicoupling_->SetScalarField(*DRT::Problem::Instance()->GetDis("structure"),scatra_->ScaTraField()->Phinp());

    // set up structural base algorithm
    struct_adapterbase_ptr_->Setup();

    // get wrapper and cast it to specific type
    // do not do so, in case the wrapper has already been set from outside
    if(structure_ == Teuchos::null)
      structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::SSIStructureWrapper>(struct_adapterbase_ptr_->StructureField());

    if(structure_ == Teuchos::null)
      dserror("No valid pointer to ADAPTER::SSIStructureWrapper !\n"
              "Either cast failed, or no valid wrapper was set using SetStructureWrapper(...) !");
  }

  // for old structural time integration
  else if(use_old_structure_)
    structure_->Setup();

  // check maps from scalar transport and structure discretizations
  if(scatra_->ScaTraField()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scalar transport discretization does not have any degrees of freedom!");
  if(structure_->DofRowMap()->NumGlobalElements() == 0)
    dserror("Structure discretization does not have any degrees of freedom!");

  // construct vector of zeroes
  zeros_ = LINALG::CreateVector(*structure_->DofRowMap());

  // set up materials
  ssicoupling_->AssignMaterialPointers(structure_->Discretization(),scatra_->ScaTraField()->Discretization());

  // set up scatra-scatra interface coupling
  if(scatra_->ScaTraField()->S2ICoupling())
  {
    // set up scatra-scatra interface coupling adapter for structure field
    SetupCouplingAdapterStructure();

    // set up map for interior and master-side structural degrees of freedom
    map_structure_condensed_ = LINALG::SplitMap(*structure_->Discretization()->DofRowMap(),*icoup_structure_->SlaveDofMap());

    // set up structural map extractor holding interior and interface maps of degrees of freedom
    std::vector<Teuchos::RCP<const Epetra_Map> > maps(0,Teuchos::null);
    maps.push_back(LINALG::SplitMap(*map_structure_condensed_,*icoup_structure_->MasterDofMap()));
    maps.push_back(icoup_structure_->SlaveDofMap());
    maps.push_back(icoup_structure_->MasterDofMap());
    maps_structure_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*structure_->Discretization()->DofRowMap(),maps));
    maps_structure_->CheckForValidMapExtractor();
  }

  // set flag
  SetIsSetup(true);

  return;
}

/*----------------------------------------------------------------------------------*
 | set up scatra-scatra interface coupling adapter for structure field   fang 08/17 |
 *----------------------------------------------------------------------------------*/
void SSI::SSI_Base::SetupCouplingAdapterStructure()
{
  // initialize integer vectors for global IDs of master-side and slave-side interface nodes on structure discretization
  std::vector<int> inodegidvec_master;
  std::vector<int> inodegidvec_slave;

  // extract scatra-scatra interface coupling conditions from structure discretization
  std::vector<DRT::Condition*> conditions(0,NULL);
  structure_->Discretization()->GetCondition("S2ICoupling",conditions);

  // loop over all conditions
  for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
  {
    // extract current condition
    DRT::Condition* const condition = conditions[icondition];

    // extract interface side associated with current condition
    const int side = condition->GetInt("interface side");

    // extract nodes associated with current condition
    const std::vector<int>* const inodegids = condition->Nodes();

    // loop over all nodes
    for(unsigned inode=0; inode<inodegids->size(); ++inode)
    {
      // extract ID of current node
      const int inodegid = (*inodegids)[inode];

      // insert global ID of current node into associated vector only if node is owned by current processor
      // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
      if(structure_->Discretization()->HaveGlobalNode(inodegid) and structure_->Discretization()->gNode(inodegid)->Owner() == structure_->Discretization()->Comm().MyPID())
        side == INPAR::S2I::side_master ? inodegidvec_master.push_back(inodegid) : inodegidvec_slave.push_back(inodegid);
    }
  }

  // remove potential duplicates from vectors
  std::sort(inodegidvec_master.begin(),inodegidvec_master.end());
  inodegidvec_master.erase(unique(inodegidvec_master.begin(),inodegidvec_master.end()),inodegidvec_master.end());
  std::sort(inodegidvec_slave.begin(),inodegidvec_slave.end());
  inodegidvec_slave.erase(unique(inodegidvec_slave.begin(),inodegidvec_slave.end()),inodegidvec_slave.end());

  // setup scatra-scatra interface coupling adapter for structure field
  icoup_structure_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoup_structure_->SetupCoupling(
      *structure_->Discretization(),
      *structure_->Discretization(),
      inodegidvec_master,
      inodegidvec_slave,
      DRT::Problem::Instance()->NDim(),
      true,
      1.e-8
      );

  return;
}

/*----------------------------------------------------------------------*
 | Setup the discretizations                                rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::InitDiscretizations(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& scatra_disname)
{
  // Scheme   : the structure discretization is received from the input. Then, an ale-scatra disc. is cloned.

  DRT::Problem* problem = DRT::Problem::Instance();

  //1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  if (scatradis->NumGlobalNodes()==0)
  {
    if(fieldcoupling_!=INPAR::SSI::coupling_volume_match)
      dserror("If 'FIELDCOUPLING' is NOT 'volume_matching' in the SSI CONTROL section cloning of the scatra discretization"
          "from the structure discretization is not supported!");

    // fill scatra discretization by cloning structure discretization
    DRT::UTILS::CloneDiscretization<SSI::ScatraStructureCloneStrategy>(structdis,scatradis);
    scatradis->FillComplete();
  }
  else
  {
    if(fieldcoupling_==INPAR::SSI::coupling_volume_match)
      dserror("Reading a TRANSPORT discretization from the .dat file for the input parameter 'FIELDCOUPLING volume_matching' in the"
          "SSI CONTROL section is not supported! As this coupling relies on matching node (and sometimes element) IDs,"
          "the ScaTra discretization is cloned from the structure discretization. Delete the ScaTra discretization"
          "from your input file.");

    // copy conditions
    // this is actually only needed for copying TRANSPORT DIRICHLET/NEUMANN CONDITIONS
    // as standard DIRICHLET/NEUMANN CONDITIONS
    std::map<std::string,std::string> conditions_to_copy;
    SSI::ScatraStructureCloneStrategy clonestrategy;
    conditions_to_copy = clonestrategy.ConditionsToCopy();
    DRT::UTILS::DiscretizationCreatorBase creator;
    creator.CopyConditions(*scatradis,*scatradis,conditions_to_copy);

    // safety check, since it is not reasonable to have SOLIDSCATRA or SOLIDPOROP1 Elements with a SCATRA::ImplType != 'impltype_undefined'
    // if they are not cloned! Therefore loop over all structure elements and check the impltype
    for(int i= 0; i< structdis->NumMyColElements(); ++i)
    {
      if(clonestrategy.GetImplType(structdis->lColElement(i)) != INPAR::SCATRA::impltype_undefined)
        dserror("A TRANSPORT discretization is read from the .dat file, which is fine since the scatra discretization is not cloned from "
            "the structure discretization. But in the STRUCTURE ELEMENTS section of the .dat file an ImplType that is NOT 'Undefined' is prescribed "
            "which does not make sense if you don't want to clone the structure discretization. Change the ImplType to 'Undefined' "
            "or decide to clone the scatra discretization from the structure discretization.");
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Setup ssi coupling object                                rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Base::InitFieldCoupling(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& scatra_disname)
{
  // initialize return variable
  int redistribution_required = (int)SSI::none;

  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  //safety check
  {
    //check for ssi coupling condition
    std::vector<DRT::Condition*> ssicoupling;
    scatradis->GetCondition("SSICoupling",ssicoupling);
    const bool havessicoupling = (ssicoupling.size()>0);

    if(havessicoupling and
        (fieldcoupling_!= INPAR::SSI::coupling_boundary_nonmatch and fieldcoupling_!= INPAR::SSI::coupling_volumeboundary_match) )
      dserror("SSICoupling condition only valid in combination with FIELDCOUPLING set to 'boundary_nonmatching' "
          "or 'volumeboundary_matching' in SSI DYNAMIC section. ");

    if(fieldcoupling_==INPAR::SSI::coupling_volume_nonmatch)
    {
      const Teuchos::ParameterList& volmortarparams = DRT::Problem::Instance()->VolmortarParams();
      if (DRT::INPUT::IntegralValue<INPAR::VOLMORTAR::CouplingType>(volmortarparams,"COUPLINGTYPE")!=
          INPAR::VOLMORTAR::couplingtype_coninter)
        dserror("Volmortar coupling only tested for consistent interpolation, "
            "i.e. 'COUPLINGTYPE consint' in VOLMORTAR COUPLING section. Try other couplings at own risk.");
    }
  }

  //build SSI coupling class
  switch(fieldcoupling_)
  {
  case INPAR::SSI::coupling_volume_match:
    ssicoupling_ = Teuchos::rcp(new SSICouplingMatchingVolume());
    break;
  case INPAR::SSI::coupling_volume_nonmatch:
    ssicoupling_ = Teuchos::rcp(new SSICouplingNonMatchingVolume());
    // redistribution is still performed inside
    redistribution_required = (int)SSI::binning;
    break;
  case INPAR::SSI::coupling_boundary_nonmatch:
    ssicoupling_ = Teuchos::rcp(new SSICouplingNonMatchingBoundary());
    break;
  case INPAR::SSI::coupling_volumeboundary_match:
    ssicoupling_ = Teuchos::rcp(new SSICouplingMatchingVolumeAndBoundary());
    redistribution_required = (int)SSI::match;
    break;
  default:
    dserror("unknown type of field coupling for SSI!");
    break;
  }

  // initialize coupling objects including dof sets
  ssicoupling_->Init(problem->NDim(),structdis,scatradis);

  return redistribution_required;
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 01/12  |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::ReadRestart( int restart )
{
  if (restart)
  {
    structure_->ReadRestart(restart);

    const Teuchos::ParameterList& ssidyn = DRT::Problem::Instance()->SSIControlParams();
    const bool restartfromstructure = DRT::INPUT::IntegralValue<int>(ssidyn,"RESTART_FROM_STRUCTURE");

    if (not restartfromstructure) //standard restart
    {
        scatra_->ScaTraField()->ReadRestart(restart);
    }
    else //restart from structure simulation
    {
      // Since there is no restart output for the scatra fiels available, we only have to fix the
      // time and step counter
      scatra_->ScaTraField()->SetTimeStep(structure_->TimeOld(),restart);
    }

    SetTimeStep(structure_->TimeOld(), restart);
  }

  // Material pointers to other field were deleted during ReadRestart().
  // They need to be reset.
  ssicoupling_->AssignMaterialPointers(structure_->Discretization(),scatra_->ScaTraField()->Discretization());

  return;
}

/*----------------------------------------------------------------------*
 | read restart information for given time (public)        AN, JH 10/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::ReadRestartfromTime( double restarttime )
{
  if ( restarttime > 0.0 )
  {
    const int restartstructure = SSI::Utils::CheckTimeStepping(structure_->Dt(), restarttime);
    const int restartscatra    = SSI::Utils::CheckTimeStepping(scatra_->ScaTraField()->Dt(), restarttime);

    structure_->ReadRestart(restartstructure);

    const Teuchos::ParameterList& ssidyn = DRT::Problem::Instance()->SSIControlParams();
    const bool restartfromstructure = DRT::INPUT::IntegralValue<int>(ssidyn,"RESTART_FROM_STRUCTURE");

    if (not restartfromstructure) //standard restart
    {
      scatra_->ScaTraField()->ReadRestart(restartscatra);
    }
    else //restart from structure simulation
    {
      // Since there is no restart output for the scatra fiels available, we only have to fix the
      // time and step counter
      scatra_->ScaTraField()->SetTimeStep(structure_->TimeOld(),restartscatra);
    }

    SetTimeStep(structure_->TimeOld(), restartstructure);
  }

  // Material pointers to other field were deleted during ReadRestart().
  // They need to be reset.
  ssicoupling_->AssignMaterialPointers(structure_->Discretization(),scatra_->ScaTraField()->Discretization());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::TestResults(const Epetra_Comm& comm) const
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
  problem->AddFieldTest(Teuchos::rcp(new SSI::SSIResultTest(Teuchos::rcp(this,false))));
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                                       Teuchos::RCP<const Epetra_Vector> vel )
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  SetMeshDisp(disp);
  SetVelocityFields(vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetScatraSolution( Teuchos::RCP<const Epetra_Vector> phi )
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  ssicoupling_->SetScalarField(*structure_->Discretization(),phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetVelocityFields( Teuchos::RCP<const Epetra_Vector> vel)
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  ssicoupling_->SetVelocityFields(scatra_,zeros_,vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp )
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  ssicoupling_->SetMeshDisp(scatra_,disp);
}


