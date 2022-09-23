/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all scalar structure algorithms

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "ssi_base.H"

#include "ssi_clonestrategy.H"
#include "ssi_coupling.H"
#include "ssi_partitioned.H"
#include "ssi_resulttest.H"
#include "ssi_str_model_evaluator_partitioned.H"
#include "ssi_utils.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"

#include "../drt_inpar/inpar_volmortar.H"
#include "../drt_inpar/inpar_ssi.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_inputreader.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/function_of_time.H"

#include "../drt_mat/matpar_bundle.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"
#include "../drt_scatra/scatra_utils.H"

#include "../drt_scatra_ele/scatra_ele.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSIBase::SSIBase(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : AlgorithmBase(comm, globaltimeparams),
      diff_time_step_size_(
          static_cast<int>(DRT::INPUT::IntegralValue<int>(globaltimeparams, "DIFFTIMESTEPSIZE"))),
      fieldcoupling_(Teuchos::getIntegralValue<INPAR::SSI::FieldCoupling>(
          DRT::Problem::Instance()->SSIControlParams(), "FIELDCOUPLING")),
      isinit_(false),
      issetup_(false),
      is_scatra_manifold_(
          DRT::INPUT::IntegralValue<bool>(globaltimeparams.sublist("MANIFOLD"), "ADD_MANIFOLD")),
      is_manifold_meshtying_(DRT::INPUT::IntegralValue<bool>(
          globaltimeparams.sublist("MANIFOLD"), "MESHTYING_MANIFOLD")),
      iter_(0),
      macro_scale_(DRT::Problem::Instance()->Materials()->FirstIdByType(
                       INPAR::MAT::m_scatra_multiscale) != -1 or
                   DRT::Problem::Instance()->Materials()->FirstIdByType(
                       INPAR::MAT::m_newman_multiscale) != -1),
      meshtying_strategy_s2i_(Teuchos::null),
      scatra_base_algorithm_(Teuchos::null),
      scatra_manifold_base_algorithm_(Teuchos::null),
      ssi_structure_meshtying_(Teuchos::null),
      ssicoupling_(Teuchos::null),
      ssiinterfacecontact_(
          DRT::Problem::Instance()->GetDis("structure")->GetCondition("SSIInterfaceContact") !=
          nullptr),
      ssiinterfacemeshtying_(
          DRT::Problem::Instance()->GetDis("structure")->GetCondition("SSIInterfaceMeshtying") !=
          nullptr),
      structure_(Teuchos::null),
      struct_adapterbase_ptr_(Teuchos::null),
      temperature_funct_num_(
          DRT::Problem::Instance()->ELCHControlParams().get<int>("TEMPERATURE_FROM_FUNCT")),
      temperature_vector_(Teuchos::null),
      use_old_structure_(false),
      zeros_structure_(Teuchos::null)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Init this class                                          rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, bool isAle)
{
  // reset the setup flag
  SetIsSetup(false);

  // do discretization specific setup (e.g. clone discr. scatra from structure)
  InitDiscretizations(comm, struct_disname, scatra_disname,
      DRT::INPUT::IntegralValue<bool>(globaltimeparams, "REDISTRIBUTE_SOLID"));

  InitTimeIntegrators(
      globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  const RedistributionType redistribution_type =
      InitFieldCoupling(comm, struct_disname, scatra_disname);

  if (redistribution_type != SSI::RedistributionType::none) Redistribute(redistribution_type);

  CheckSSIFlags();

  CheckSSIInterfaceConditions(struct_disname);

  // set isinit_ flag true
  SetIsInit(true);
}

/*----------------------------------------------------------------------*
 | Setup this class                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::Setup()
{
  // check initialization
  CheckIsInit();

  // set up helper class for field coupling
  ssicoupling_->Setup();

  // in case of an ssi  multi scale formulation we need to set the displacement here
  auto dummy_vec = Teuchos::rcp(
      new Epetra_Vector(*DRT::Problem::Instance()->GetDis("structure")->DofRowMap(), true));
  ssicoupling_->SetMeshDisp(ScaTraBaseAlgorithm(), dummy_vec);

  // set up scalar transport field
  ScaTraField()->Setup();
  if (IsScaTraManifold()) ScaTraManifold()->Setup();

  // only relevant for new structural time integration
  // only if adapter base has not already been set up outside
  if (not use_old_structure_ and not struct_adapterbase_ptr_->IsSetup())
  {
    // set up structural model evaluator
    SetupModelEvaluator();

    // pass initial scalar field to structural discretization to correctly compute initial
    // accelerations
    if (Teuchos::getIntegralValue<INPAR::SSI::SolutionSchemeOverFields>(
            DRT::Problem::Instance()->SSIControlParams(), "COUPALGO") !=
        INPAR::SSI::SolutionSchemeOverFields::ssi_OneWay_SolidToScatra)
      ssicoupling_->SetScalarField(
          *DRT::Problem::Instance()->GetDis("structure"), ScaTraField()->Phinp(), 1);

    if (macro_scale_)
    {
      ScaTraField()->CalcMeanMicroConcentration();
      ssicoupling_->SetScalarFieldMicro(
          *DRT::Problem::Instance()->GetDis("structure"), ScaTraField()->PhinpMicro(), 2);
    }


    //   temperature is non primary variable. Only set, if function for temperature is given
    if (temperature_funct_num_ != -1)
    {
      temperature_vector_ = Teuchos::rcp(
          new Epetra_Vector(*DRT::Problem::Instance()->GetDis("structure")->DofRowMap(2), true));

      temperature_vector_->PutScalar(
          DRT::Problem::Instance()
              ->FunctionById<DRT::UTILS::FunctionOfTime>(temperature_funct_num_ - 1)
              .Evaluate(Time()));

      ssicoupling_->SetTemperatureField(
          *DRT::Problem::Instance()->GetDis("structure"), temperature_vector_);
    }

    // set up structural base algorithm
    struct_adapterbase_ptr_->Setup();

    // get wrapper and cast it to specific type
    // do not do so, in case the wrapper has already been set from outside
    if (structure_ == Teuchos::null)
      structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::SSIStructureWrapper>(
          struct_adapterbase_ptr_->StructureField());

    if (structure_ == Teuchos::null)
    {
      dserror(
          "No valid pointer to ADAPTER::SSIStructureWrapper !\n"
          "Either cast failed, or no valid wrapper was set using SetStructureWrapper(...) !");
    }
  }

  // for old structural time integration
  else if (use_old_structure_)
    structure_->Setup();

  // check maps from scalar transport and structure discretizations
  if (ScaTraField()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scalar transport discretization does not have any degrees of freedom!");
  if (structure_->DofRowMap()->NumGlobalElements() == 0)
    dserror("Structure discretization does not have any degrees of freedom!");

  // set up materials
  ssicoupling_->AssignMaterialPointers(
      structure_->Discretization(), ScaTraField()->Discretization());

  // set up scatra-scatra interface coupling
  if (SSIInterfaceMeshtying())
  {
    ssi_structure_meshtying_ = Teuchos::rcp(new SSI::UTILS::SSIStructureMeshTying(
        "SSIInterfaceMeshtying", structure_->Discretization()));

    // extract meshtying strategy for scatra-scatra interface coupling on scatra discretization
    meshtying_strategy_s2i_ =
        Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(ScaTraField()->Strategy());

    // safety checks
    if (meshtying_strategy_s2i_ == Teuchos::null)
      dserror("Invalid scatra-scatra interface coupling strategy!");
  }

  // construct vector of zeroes
  zeros_structure_ = LINALG::CreateVector(*structure_->DofRowMap());

  // set flag
  SetIsSetup(true);
}

/*----------------------------------------------------------------------*
 | Setup the discretizations                                rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::InitDiscretizations(const Epetra_Comm& comm, const std::string& struct_disname,
    const std::string& scatra_disname, bool redistribute_struct_dis)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  auto structdis = problem->GetDis(struct_disname);
  auto scatradis = problem->GetDis(scatra_disname);

  if (redistribute_struct_dis) DRT::UTILS::RedistributeDiscretizationsByBinning({structdis}, false);

  if (scatradis->NumGlobalNodes() == 0)
  {
    if (fieldcoupling_ != INPAR::SSI::FieldCoupling::volume_match and
        fieldcoupling_ != INPAR::SSI::FieldCoupling::volumeboundary_match)
    {
      dserror(
          "If 'FIELDCOUPLING' is NOT 'volume_matching' or 'volumeboundary_matching' in the SSI "
          "CONTROL section cloning of the scatra discretization from the structure discretization "
          "is not supported!");
    }

    // fill scatra discretization by cloning structure discretization
    DRT::UTILS::CloneDiscretization<SSI::ScatraStructureCloneStrategy>(structdis, scatradis);
    scatradis->FillComplete();

    // create discretization for scatra manifold based on SSISurfaceManifold condition
    if (IsScaTraManifold())
    {
      auto scatra_manifold_dis = problem->GetDis("scatra_manifold");
      DRT::UTILS::CloneDiscretizationFromCondition<SSI::ScatraStructureCloneStrategyManifold>(
          *structdis, *scatra_manifold_dis, "SSISurfaceManifold");

      // clone conditions. Needed this way, as many conditions are cloned from SSISurfaceManifold.
      std::vector<std::map<std::string, std::string>> conditions_to_copy = {
          {std::make_pair("SSISurfaceManifold", "SSISurfaceManifold")},
          {std::make_pair("ScaTraManifoldInitfield", "Initfield")},
          {std::make_pair("ManifoldDirichlet", "Dirichlet")}};

      //! in case of no mesh tying between manifolds: partition manifold domains
      if (!is_manifold_meshtying_)
      {
        std::map<std::string, std::string> temp_map = {
            std::make_pair("SSISurfaceManifold", "ScatraPartitioning")};
        conditions_to_copy.emplace_back(temp_map);
      }

      const auto output_scalar_type = DRT::INPUT::IntegralValue<INPAR::SCATRA::OutputScalarType>(
          problem->ScalarTransportDynamicParams(), "OUTPUTSCALARS");
      if (output_scalar_type == INPAR::SCATRA::outputscalars_condition or
          output_scalar_type == INPAR::SCATRA::outputscalars_entiredomain_condition)
      {
        std::map<std::string, std::string> tempmap = {
            std::make_pair("SSISurfaceManifold", "TotalAndMeanScalar")};

        conditions_to_copy.emplace_back(tempmap);
      }

      DRT::UTILS::DiscretizationCreatorBase creator;
      for (const auto& condition_to_copy : conditions_to_copy)
        creator.CopyConditions(*structdis, *scatra_manifold_dis, condition_to_copy);

      scatra_manifold_dis->FillComplete();

      //! in case of mesh tying between manifolds: unite manifold domains -> create new
      //! ScatraPartitioning condition
      if (is_manifold_meshtying_)
      {
        // create vector of all node GIDs (all procs) of manifold dis
        int num_my_nodes = scatra_manifold_dis->NodeRowMap()->NumMyElements();
        std::vector<int> my_node_ids(num_my_nodes);
        for (int lid = 0; lid < num_my_nodes; ++lid)
          my_node_ids[lid] = scatra_manifold_dis->NodeRowMap()->GID(lid);

        int max_num_nodes = 0;
        Comm().MaxAll(&num_my_nodes, &max_num_nodes, 1);

        // resize vector and fill with place holders (-1)
        my_node_ids.resize(max_num_nodes, -1);

        std::vector<int> glob_node_ids(max_num_nodes * Comm().NumProc(), -1);
        Comm().GatherAll(&my_node_ids[0], &glob_node_ids[0], static_cast<int>(my_node_ids.size()));

        // remove place holders (-1)
        glob_node_ids.erase(
            std::remove(glob_node_ids.begin(), glob_node_ids.end(), -1), glob_node_ids.end());

        // create new condition
        const int num_conditions = static_cast<int>(scatra_manifold_dis->GetAllConditions().size());
        auto cond = Teuchos::rcp(new DRT::Condition(
            num_conditions + 1, DRT::Condition::ScatraPartitioning, true, DRT::Condition::Surface));
        cond->Add("ConditionID", 0);
        cond->Add("Node Ids", glob_node_ids);

        scatra_manifold_dis->SetCondition("ScatraPartitioning", cond);

        scatra_manifold_dis->FillComplete();
      }
    }
  }
  else
  {
    if (fieldcoupling_ == INPAR::SSI::FieldCoupling::volume_match)
    {
      dserror(
          "Reading a TRANSPORT discretization from the .dat file for the input parameter "
          "'FIELDCOUPLING volume_matching' in the SSI CONTROL section is not supported! As this "
          "coupling relies on matching node (and sometimes element) IDs, the ScaTra discretization "
          "is cloned from the structure discretization. Delete the ScaTra discretization from your "
          "input file.");
    }

    // copy conditions
    // this is actually only needed for copying TRANSPORT DIRICHLET/NEUMANN CONDITIONS
    // as standard DIRICHLET/NEUMANN CONDITIONS
    std::map<std::string, std::string> conditions_to_copy;
    SSI::ScatraStructureCloneStrategy clonestrategy;
    conditions_to_copy = clonestrategy.ConditionsToCopy();
    DRT::UTILS::DiscretizationCreatorBase creator;
    creator.CopyConditions(*scatradis, *scatradis, conditions_to_copy);

    // safety check, since it is not reasonable to have SOLIDSCATRA or SOLIDPOROP1 Elements with a
    // SCATRA::ImplType != 'impltype_undefined' if they are not cloned! Therefore loop over all
    // structure elements and check the impltype
    for (int i = 0; i < structdis->NumMyColElements(); ++i)
    {
      if (clonestrategy.GetImplType(structdis->lColElement(i)) != INPAR::SCATRA::impltype_undefined)
      {
        dserror(
            "A TRANSPORT discretization is read from the .dat file, which is fine since the scatra "
            "discretization is not cloned from the structure discretization. But in the STRUCTURE "
            "ELEMENTS section of the .dat file an ImplType that is NOT 'Undefined' is prescribed "
            "which does not make sense if you don't want to clone the structure discretization. "
            "Change the ImplType to 'Undefined' or decide to clone the scatra discretization from "
            "the structure discretization.");
      }
    }
  }
  // read in the micro field, has to be done after cloning of the scatra discretization
  auto input_file_name = problem->OutputControlFile()->InputFileName();
  DRT::INPUT::DatFileReader local_reader(input_file_name);
  problem->ReadMicroFields(local_reader);
}

/*----------------------------------------------------------------------*
 | Setup ssi coupling object                                rauch 08/16 |
 *----------------------------------------------------------------------*/
SSI::RedistributionType SSI::SSIBase::InitFieldCoupling(
    const Epetra_Comm& comm, const std::string& struct_disname, const std::string& scatra_disname)
{
  // initialize return variable
  RedistributionType redistribution_required = SSI::RedistributionType::none;

  DRT::Problem* problem = DRT::Problem::Instance();
  auto structdis = problem->GetDis(struct_disname);
  auto scatradis = problem->GetDis(scatra_disname);
  Teuchos::RCP<DRT::Discretization> scatra_manifold_dis(Teuchos::null);
  if (IsScaTraManifold()) scatra_manifold_dis = problem->GetDis("scatra_manifold");

  // safety check
  {
    // check for ssi coupling condition
    std::vector<DRT::Condition*> ssicoupling;
    scatradis->GetCondition("SSICoupling", ssicoupling);
    const bool havessicoupling = (ssicoupling.size() > 0);

    if (havessicoupling and (fieldcoupling_ != INPAR::SSI::FieldCoupling::boundary_nonmatch and
                                fieldcoupling_ != INPAR::SSI::FieldCoupling::volumeboundary_match))
    {
      dserror(
          "SSICoupling condition only valid in combination with FIELDCOUPLING set to "
          "'boundary_nonmatching' or 'volumeboundary_matching' in SSI DYNAMIC section. ");
    }

    if (fieldcoupling_ == INPAR::SSI::FieldCoupling::volume_nonmatch)
    {
      const Teuchos::ParameterList& volmortarparams = DRT::Problem::Instance()->VolmortarParams();
      if (DRT::INPUT::IntegralValue<INPAR::VOLMORTAR::CouplingType>(
              volmortarparams, "COUPLINGTYPE") != INPAR::VOLMORTAR::couplingtype_coninter)
      {
        dserror(
            "Volmortar coupling only tested for consistent interpolation, "
            "i.e. 'COUPLINGTYPE consint' in VOLMORTAR COUPLING section. Try other couplings at own "
            "risk.");
      }
    }
    if (IsScaTraManifold() and fieldcoupling_ != INPAR::SSI::FieldCoupling::volumeboundary_match)
      dserror("Solving manifolds only in combination with matching volumes and boundaries");
  }

  // build SSI coupling class
  switch (fieldcoupling_)
  {
    case INPAR::SSI::FieldCoupling::volume_match:
      ssicoupling_ = Teuchos::rcp(new SSICouplingMatchingVolume());
      break;
    case INPAR::SSI::FieldCoupling::volume_nonmatch:
      ssicoupling_ = Teuchos::rcp(new SSICouplingNonMatchingVolume());
      // redistribution is still performed inside
      redistribution_required = SSI::RedistributionType::binning;
      break;
    case INPAR::SSI::FieldCoupling::boundary_nonmatch:
      ssicoupling_ = Teuchos::rcp(new SSICouplingNonMatchingBoundary());
      break;
    case INPAR::SSI::FieldCoupling::volumeboundary_match:
      ssicoupling_ = Teuchos::rcp(new SSICouplingMatchingVolumeAndBoundary());
      redistribution_required = SSI::RedistributionType::match;
      break;
    default:
      dserror("unknown type of field coupling for SSI!");
      break;
  }

  // initialize coupling objects including dof sets
  ssicoupling_->Init(problem->NDim(), structdis, scatradis, scatra_manifold_dis);

  return redistribution_required;
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 01/12  |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::ReadRestart(int restart)
{
  if (restart)
  {
    structure_->ReadRestart(restart);

    const Teuchos::ParameterList& ssidyn = DRT::Problem::Instance()->SSIControlParams();
    const bool restartfromstructure =
        DRT::INPUT::IntegralValue<int>(ssidyn, "RESTART_FROM_STRUCTURE");

    if (not restartfromstructure)  // standard restart
    {
      ScaTraField()->ReadRestart(restart);
      if (IsScaTraManifold()) ScaTraManifold()->ReadRestart(restart);
    }
    else  // restart from structure simulation
    {
      // Since there is no restart output for the scatra field available, we only have to fix the
      // time and step counter
      ScaTraField()->SetTimeStep(structure_->TimeOld(), restart);
      if (IsScaTraManifold()) ScaTraManifold()->SetTimeStep(structure_->TimeOld(), restart);
    }

    SetTimeStep(structure_->TimeOld(), restart);
  }

  // Material pointers to other field were deleted during ReadRestart().
  // They need to be reset.
  ssicoupling_->AssignMaterialPointers(
      structure_->Discretization(), ScaTraField()->Discretization());
}

/*----------------------------------------------------------------------*
 | read restart information for given time (public)        AN, JH 10/14 |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::ReadRestartfromTime(double restarttime)
{
  if (restarttime > 0.0)
  {
    const int restartstructure = SSI::UTILS::CheckTimeStepping(structure_->Dt(), restarttime);
    const int restartscatra = SSI::UTILS::CheckTimeStepping(ScaTraField()->Dt(), restarttime);
    const int restartscatramanifold =
        IsScaTraManifold() ? SSI::UTILS::CheckTimeStepping(ScaTraManifold()->Dt(), restarttime)
                           : -1;

    structure_->ReadRestart(restartstructure);

    const Teuchos::ParameterList& ssidyn = DRT::Problem::Instance()->SSIControlParams();
    const bool restartfromstructure =
        DRT::INPUT::IntegralValue<int>(ssidyn, "RESTART_FROM_STRUCTURE");

    if (not restartfromstructure)  // standard restart
    {
      ScaTraField()->ReadRestart(restartscatra);
      if (IsScaTraManifold()) ScaTraManifold()->ReadRestart(restartscatramanifold);
    }
    else  // restart from structure simulation
    {
      // Since there is no restart output for the scatra fiels available, we only have to fix the
      // time and step counter
      ScaTraField()->SetTimeStep(structure_->TimeOld(), restartscatra);
      if (IsScaTraManifold())
        ScaTraManifold()->SetTimeStep(structure_->TimeOld(), restartscatramanifold);
    }

    SetTimeStep(structure_->TimeOld(), restartstructure);
  }

  // Material pointers to other field were deleted during ReadRestart().
  // They need to be reset.
  ssicoupling_->AssignMaterialPointers(
      structure_->Discretization(), ScaTraField()->Discretization());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::TestResults(const Epetra_Comm& comm) const
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(ScaTraBaseAlgorithm()->CreateScaTraFieldTest());
  if (IsScaTraManifold())
    problem->AddFieldTest(ScaTraManifoldBaseAlgorithm()->CreateScaTraFieldTest());
  problem->AddFieldTest(Teuchos::rcp(new SSI::SSIResultTest(Teuchos::rcp(this, false))));
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetStructSolution(
    Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel)
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  SetMeshDisp(disp);
  SetVelocityFields(vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  ssicoupling_->SetScalarField(*StructureField()->Discretization(), phi, 1);
  if (IsScaTraManifold()) ssicoupling_->SetScalarField(*ScaTraManifold()->Discretization(), phi, 2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetMicroScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  ssicoupling_->SetScalarFieldMicro(*StructureField()->Discretization(), phi, 2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetScatraManifoldSolution(Teuchos::RCP<const Epetra_Vector> phi)
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();
  ssicoupling_->SetScaTraManifoldField(*ScaTraField()->Discretization(), phi, 2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::EvaluateAndSetTemperatureField()
{
  // temperature is non primary variable. Only set, if function for temperature is given
  if (temperature_funct_num_ != -1)
  {
    // evaluate temperature at current time and put to scalar
    const double temperature =
        DRT::Problem::Instance()
            ->FunctionById<DRT::UTILS::FunctionOfTime>(temperature_funct_num_ - 1)
            .Evaluate(Time());
    temperature_vector_->PutScalar(temperature);

    // set temperature vector to structure discretization
    ssicoupling_->SetTemperatureField(*structure_->Discretization(), temperature_vector_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetVelocityFields(Teuchos::RCP<const Epetra_Vector> vel)
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  ssicoupling_->SetVelocityFields(ScaTraBaseAlgorithm(), zeros_structure_, vel);
  if (IsScaTraManifold())
    ssicoupling_->SetVelocityFields(ScaTraManifoldBaseAlgorithm(), zeros_structure_, vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetMeshDisp(Teuchos::RCP<const Epetra_Vector> disp)
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  ssicoupling_->SetMeshDisp(ScaTraBaseAlgorithm(), disp);
  if (IsScaTraManifold()) ssicoupling_->SetMeshDisp(ScaTraManifoldBaseAlgorithm(), disp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::CheckSSIFlags() const
{
  if (ScaTraField()->S2IKinetics())
  {
    if (!(SSIInterfaceContact() or SSIInterfaceMeshtying()))
    {
      dserror(
          "You defined an 'S2IKinetics' condition in the input-file. However, neither an "
          "'SSIInterfaceContact' condition nor an 'SSIInterfaceMeshtying' condition defined. This "
          "is not reasonable!");
    }
  }

  const bool is_nitsche_penalty_adaptive(DRT::INPUT::IntegralValue<int>(
      DRT::Problem::Instance()->ContactDynamicParams(), "NITSCHE_PENALTY_ADAPTIVE"));

  if (SSIInterfaceContact() and is_nitsche_penalty_adaptive)
    dserror("Adaptive nitsche penalty parameter currently not supported!");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetDtFromScaTraToStructure()
{
  StructureField()->SetDt(ScaTraField()->Dt());
  StructureField()->SetTimen(ScaTraField()->Time());
  StructureField()->PostUpdate();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetDtFromScaTraToManifold()
{
  ScaTraManifold()->SetDt(ScaTraField()->Dt());
  ScaTraManifold()->SetTimeStep(ScaTraField()->Time(), ScaTraField()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetDtFromScaTraToSSI()
{
  // set values for this SSI algorithm
  SetTimeStep(ScaTraField()->Time(), Step());
  SetDt(ScaTraField()->Dt());

  // set values for other fields
  SetDtFromScaTraToStructure();
  if (IsScaTraManifold()) SetDtFromScaTraToManifold();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::Redistribute(const RedistributionType redistribution_type)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  auto structdis = problem->GetDis("structure");
  auto scatradis = problem->GetDis("scatra");
  if (redistribution_type == SSI::RedistributionType::match and !IsScaTraManifold())
  {
    // first we bin the scatra discretization
    std::vector<Teuchos::RCP<DRT::Discretization>> dis;
    dis.push_back(scatradis);
    DRT::UTILS::RedistributeDiscretizationsByBinning(dis, false);

    DRT::UTILS::MatchElementDistributionOfMatchingConditionedElements(
        *scatradis, *scatradis, "ScatraHeteroReactionMaster", "ScatraHeteroReactionSlave");

    // now we redistribute the structure dis to match the scatra dis
    DRT::UTILS::MatchElementDistributionOfMatchingDiscretizations(*scatradis, *structdis);
  }
  else if (redistribution_type == SSI::RedistributionType::binning)
  {
    // create vector of discr.
    std::vector<Teuchos::RCP<DRT::Discretization>> dis;
    dis.push_back(structdis);
    dis.push_back(scatradis);

    DRT::UTILS::RedistributeDiscretizationsByBinning(dis, false);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> SSI::SSIBase::ScaTraField() const
{
  return scatra_base_algorithm_->ScaTraField();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> SSI::SSIBase::ScaTraManifold() const
{
  return scatra_manifold_base_algorithm_->ScaTraField();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::InitTimeIntegrators(const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, const bool isAle)
{
  // get the global problem
  auto* problem = DRT::Problem::Instance();

  // time parameter handling
  // In case of different time stepping, time params have to be read from single field sections.
  // In case of equal time step size for all fields the time params are controlled solely by the
  // problem section (e.g. ---SSI DYNAMIC or ---CELL DYNAMIC).
  const auto* structtimeparams = &globaltimeparams;
  const auto* scatratimeparams = &globaltimeparams;
  if (diff_time_step_size_)
  {
    structtimeparams = &structparams;
    scatratimeparams = &scatraparams;
  }

  // we do not construct a structure, in case it was built externally and handed into this object
  if (struct_adapterbase_ptr_ == Teuchos::null)
  {
    // access the structural discretization
    auto structdis = problem->GetDis(struct_disname);

    // build structure based on new structural time integration
    if (structparams.get<std::string>("INT_STRATEGY") == "Standard")
    {
      struct_adapterbase_ptr_ = ADAPTER::STR::BuildStructureAlgorithm(structparams);

      // initialize structure base algorithm
      struct_adapterbase_ptr_->Init(
          *structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis);
    }
    // build structure based on old structural time integration
    else if (structparams.get<std::string>("INT_STRATEGY") == "Old")
    {
      use_old_structure_ = true;

      auto structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
          *structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis));
      structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::SSIStructureWrapper>(
          structure->StructureField(), true);
      if (structure_ == Teuchos::null)
        dserror("cast from ADAPTER::Structure to ADAPTER::SSIStructureWrapper failed");
    }
    else
    {
      dserror(
          "Unknown time integration requested!\n"
          "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!\n"
          "If you want to use yet unsupported elements or algorithms,\n"
          "set INT_STRATEGY to Old in ---STRUCUTRAL DYNAMIC section!");
    }
  }

  // create and initialize scatra base algorithm.
  // scatra time integrator constructed and initialized inside.
  // mesh is written inside. cloning must happen before!
  scatra_base_algorithm_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

  ScaTraBaseAlgorithm()->Init(*scatratimeparams, SSI::UTILS::ModifyScaTraParams(scatraparams),
      problem->SolverParams(scatraparams.get<int>("LINEAR_SOLVER")), scatra_disname, isAle);
  ScaTraBaseAlgorithm()->ScaTraField()->SetNumberOfDofSetDisplacement(1);
  ScaTraBaseAlgorithm()->ScaTraField()->SetNumberOfDofSetVelocity(1);

  if (macro_scale_) ScaTraBaseAlgorithm()->ScaTraField()->SetNumberOfDofSetMicroScale(2);

  // create and initialize scatra base algorithm for manifolds
  if (IsScaTraManifold())
  {
    scatra_manifold_base_algorithm_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

    ScaTraManifoldBaseAlgorithm()->Init(*scatratimeparams,
        SSI::UTILS::CloneScaTraManifoldParams(scatraparams, globaltimeparams.sublist("MANIFOLD")),
        problem->SolverParams(globaltimeparams.sublist("MANIFOLD").get<int>("LINEAR_SOLVER")),
        "scatra_manifold", isAle);
    ScaTraManifoldBaseAlgorithm()->ScaTraField()->SetNumberOfDofSetDisplacement(1);
    ScaTraManifoldBaseAlgorithm()->ScaTraField()->SetNumberOfDofSetVelocity(1);
  }

  // do checks if adaptive time stepping is activated
  if (DRT::INPUT::IntegralValue<bool>(globaltimeparams, "ADAPTIVE_TIMESTEPPING"))
    CheckAdaptiveTimeStepping(scatraparams, structparams);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SSI::SSIBase::DoCalculateInitialPotentialField() const
{
  const auto ssi_params = DRT::Problem::Instance()->SSIControlParams();
  const bool init_pot_calc =
      DRT::INPUT::IntegralValue<bool>(ssi_params.sublist("ELCH"), "INITPOTCALC");

  return init_pot_calc and IsElchScaTraTimIntType();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SSI::SSIBase::IsElchScaTraTimIntType() const
{
  const auto ssi_params = DRT::Problem::Instance()->SSIControlParams();
  const auto scatra_type =
      Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(ssi_params, "SCATRATIMINTTYPE");

  return scatra_type == INPAR::SSI::ScaTraTimIntType::elch;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SSI::SSIBase::IsRestart() const
{
  // get the global problem
  const auto* problem = DRT::Problem::Instance();

  const int restartstep = problem->Restart();
  const double restarttime = problem->RestartTime();

  return ((restartstep > 0) or (restarttime > 0.0));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::CheckAdaptiveTimeStepping(
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams)
{
  // safety check: adaptive time stepping in one of the sub problems
  if (!DRT::INPUT::IntegralValue<bool>(scatraparams, "ADAPTIVE_TIMESTEPPING"))
  {
    dserror(
        "Must provide adaptive time stepping algorithm in one of the sub problems. (Currently "
        "just ScaTra)");
  }
  if (DRT::INPUT::IntegralValue<int>(structparams.sublist("TIMEADAPTIVITY"), "KIND") !=
      INPAR::STR::timada_kind_none)
    dserror("Adaptive time stepping in SSI currently just from ScaTra");
  if (DRT::INPUT::IntegralValue<int>(structparams, "DYNAMICTYP") == INPAR::STR::dyna_ab2)
    dserror("Currently, only one step methods are allowed for adaptive time stepping");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::CheckSSIInterfaceConditions(const std::string& struct_disname) const
{
  // access the structural discretization
  auto structdis = DRT::Problem::Instance()->GetDis(struct_disname);

  if (SSIInterfaceMeshtying())
    SCATRA::SCATRAUTILS::CheckConsistencyWithS2IKineticsCondition(
        "SSIInterfaceMeshtying", structdis);

  // check scatra-structure-interaction contact condition
  if (SSIInterfaceContact())
  {
    // get ssi condition to be tested
    std::vector<DRT::Condition*> ssiconditions;
    structdis->GetCondition("SSIInterfaceContact", ssiconditions);
    SSI::UTILS::CheckConsistencyOfSSIInterfaceContactCondition(ssiconditions, structdis);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetupSystem()
{
  if (ssiinterfacemeshtying_)
    SSIStructureMeshTying()->CheckSlaveSideHasDirichletConditions(
        StructureField()->GetDBCMapExtractor()->CondMap());
}