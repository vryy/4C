/*----------------------------------------------------------------------*/
/*! \file
\brief General algorithmic routines for partitioned solution approaches
       to fluid-structure-scalar-scalar interaction (FS3I), that is,
       algorithmic routines not specifically related to partitioned
       solution approaches to one -or two-way-coupled problem
       configurations, respectively

\level 2



*----------------------------------------------------------------------*/


#include "4C_fs3i_partitioned.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_structure_scatra_ele.hpp"
#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fem_condition_selector.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_monolithicfluidsplit.hpp"
#include "4C_fsi_monolithicstructuresplit.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"
#include "4C_ssi_clonestrategy.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I::PartFS3I(const Epetra_Comm& comm) : FS3IBase(), comm_(comm)
{
  // Keep constructor empty!
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::init()
{
  // call setup in base class
  FS3I::FS3IBase::init();

  volume_fieldcouplings_.push_back(Core::UTILS::IntegralValue<Inpar::FS3I::VolumeCoupling>(
      Global::Problem::Instance()->FS3IDynamicParams(), "FLUIDSCAL_FIELDCOUPLING"));
  volume_fieldcouplings_.push_back(Core::UTILS::IntegralValue<Inpar::FS3I::VolumeCoupling>(
      Global::Problem::Instance()->FS3IDynamicParams(), "STRUCTSCAL_FIELDCOUPLING"));

  Global::Problem* problem = Global::Problem::Instance();
  const Teuchos::ParameterList& fs3idyn = problem->FS3IDynamicParams();

  //---------------------------------------------------------------------
  // ensure correct order of three discretizations, with dof-numbering
  // such that structure dof < fluid dof < ale dofs
  // (ordering required at certain non-intuitive points)
  //---------------------------------------------------------------------
  problem->GetDis("structure")->fill_complete();
  problem->GetDis("fluid")->fill_complete();
  problem->GetDis("ale")->fill_complete();
  problem->GetDis("scatra1")->fill_complete();
  problem->GetDis("scatra2")->fill_complete();

  //---------------------------------------------------------------------
  // access discretizations for structure, fluid, ale as well as fluid-
  // and structure-based scalar transport and get material map for fluid
  // and scalar transport elements
  //---------------------------------------------------------------------
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<Core::FE::Discretization> fluidscatradis = problem->GetDis("scatra1");
  Teuchos::RCP<Core::FE::Discretization> structscatradis = problem->GetDis("scatra2");
  Teuchos::RCP<Core::FE::Discretization> aledis = problem->GetDis("ale");

  //---------------------------------------------------------------------
  // create ale discretization as a clone from fluid discretization
  //---------------------------------------------------------------------
  if (aledis->NumGlobalNodes() == 0)
  {
    Core::FE::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(
        fluiddis, aledis, Global::Problem::Instance()->CloningMaterialMap());
    aledis->fill_complete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->evaluate(params);
  }
  else
    FOUR_C_THROW("Providing an ALE mesh is not supported for FS3I problems.");

  // std::map<std::pair<std::string,std::string>,std::map<int,int> > clonefieldmatmap =
  // problem->CloningMaterialMap(); if (clonefieldmatmap.size() < 2)
  //  FOUR_C_THROW("At least two material lists required for partitioned FS3I!");

  // determine type of scalar transport
  const Inpar::ScaTra::ImplType impltype_fluid =
      Core::UTILS::IntegralValue<Inpar::ScaTra::ImplType>(
          Global::Problem::Instance()->FS3IDynamicParams(), "FLUIDSCAL_SCATRATYPE");

  //---------------------------------------------------------------------
  // create discretization for fluid-based scalar transport from and
  // according to fluid discretization
  //---------------------------------------------------------------------
  if (fluiddis->NumGlobalNodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

  // create fluid-based scalar transport elements if fluid-based scalar
  // transport discretization is empty
  if (fluidscatradis->NumGlobalNodes() == 0)
  {
    if (not(volume_fieldcouplings_[0] == Inpar::FS3I::coupling_match))
      FOUR_C_THROW(
          "If you clone your fluid-scatra mesh from the fluid use FLUIDSCAL_FIELDCOUPLING "
          "'volume_matching'!");

    // fill fluid-based scatra discretization by cloning fluid discretization
    Core::FE::CloneDiscretization<ScaTra::ScatraFluidCloneStrategy>(
        fluiddis, fluidscatradis, Global::Problem::Instance()->CloningMaterialMap());
    fluidscatradis->fill_complete();
    // set implementation type of cloned scatra elements to advanced reactions
    for (int i = 0; i < fluidscatradis->NumMyColElements(); ++i)
    {
      Discret::ELEMENTS::Transport* element =
          dynamic_cast<Discret::ELEMENTS::Transport*>(fluidscatradis->lColElement(i));
      if (element == nullptr)
        FOUR_C_THROW("Invalid element type!");
      else
        element->SetImplType(impltype_fluid);
    }

    volume_coupling_objects_.push_back(Teuchos::null);

    // care for secondary dof sets:
    // add proxy of fluid degrees of freedom to scatra discretization
    if (fluidscatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
      FOUR_C_THROW("Fluid scatra discretization has illegal number of dofsets!");
  }
  else
  {
    if (not(volume_fieldcouplings_[0] == Inpar::FS3I::coupling_nonmatch))
      FOUR_C_THROW(
          "If you have specified the fluid-scalar by TRANSPORT ELEMENTS use "
          "FLUIDSCAL_FIELDCOUPLING 'volume_nonmatching'!");

    if (not(impltype_fluid == Inpar::ScaTra::impltype_undefined))
      FOUR_C_THROW(
          "Be aware that your FLUIDSCAL_SCATRATYPE will be ignored and the impltype from the "
          "TRANSPORT ELMENTS section will be utilized. Use FLUIDSCAL_SCATRATYPE 'Undefined'!");

    volume_coupling_objects_.push_back(create_vol_mortar_object(fluiddis, fluidscatradis));

    FOUR_C_THROW("Mortar volume coupling for the fluid-scalar is yet not tested. So be careful!");
  }

  //---------------------------------------------------------------------
  // create discretization for structure-based scalar transport from and
  // according to structure discretization
  //---------------------------------------------------------------------
  if (structdis->NumGlobalNodes() == 0) FOUR_C_THROW("Structure discretization is empty!");

  // create structure-based scalar transport elements if structure-based
  // scalar transport discretization is empty
  if (structscatradis->NumGlobalNodes() == 0)
  {
    if (not(volume_fieldcouplings_[1] == Inpar::FS3I::coupling_match))
      FOUR_C_THROW(
          "If you clone your structure-scatra mesh from the structure use STRUCTSCAL_FIELDCOUPLING "
          "'volume_matching'!");

    // fill structure-based scatra discretization by cloning structure discretization
    Core::FE::CloneDiscretization<SSI::ScatraStructureCloneStrategy>(
        structdis, structscatradis, Global::Problem::Instance()->CloningMaterialMap());
    structscatradis->fill_complete();

    volume_coupling_objects_.push_back(Teuchos::null);

    // care for secondary dof sets:
    // add proxy of structure scatra degrees of freedom to structure discretization
    if (structdis->AddDofSet(structscatradis->GetDofSetProxy()) != 1)
      FOUR_C_THROW("Structure discretization has illegal number of dofsets!");

    // add proxy of structure degrees of freedom to scatra discretization
    if (structscatradis->AddDofSet(structdis->GetDofSetProxy()) != 1)
      FOUR_C_THROW("Structure scatra discretization has illegal number of dofsets!");
  }
  else
  {
    if (not(volume_fieldcouplings_[1] == Inpar::FS3I::coupling_nonmatch))
      FOUR_C_THROW(
          "If you have specified the structure-scalar by TRANSPORT2 ELEMENTS use "
          "STRUCTSCAL_FIELDCOUPLING 'volume_nonmatching'!");

    // is the set ImplType for the STRUCTURE Elements reasonable in case they are not cloned?
    for (int i = 0; i < structdis->NumMyColElements(); ++i)
    {
      if (Adapter::GetScaTraImplType(structdis->lColElement(i)) !=
          Inpar::ScaTra::impltype_undefined)
      {
        FOUR_C_THROW(
            "Be aware that the ImplType defined for the STRUCTURE Elements will be ignored and the "
            "ImplType from the TRANSPORT2 ELMENTS section will be utilized. Use TYPE 'Undefined' "
            "if "
            "cloning the scatra discretization from structure discretization is not intended!");
      }
    }

    volume_coupling_objects_.push_back(create_vol_mortar_object(structdis, structscatradis));
  }

  // safety check
  if (not(volume_coupling_objects_.size() == 2))
    FOUR_C_THROW("Unexpected size of volmortar object vector!");

  fluiddis->fill_complete(true, false, false);
  structdis->fill_complete(true, false, false);
  fluidscatradis->fill_complete(true, false, false);
  structscatradis->fill_complete(true, false, false);

  // Note: in the scatra fields we have now the following dof-sets:
  // structure dofset 0: structure dofset
  // structure dofset 1: structscatra dofset
  // fluidscatra dofset 0: fluidscatra dofset
  // fluidscatra dofset 1: fluid dofset
  // structscatra dofset 0: structscatra dofset
  // structscatra dofset 1: structure dofset

  //---------------------------------------------------------------------
  // get FSI coupling algorithm
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  int coupling = Teuchos::getIntegralValue<int>(fsidyn, "COUPALGO");

  const Teuchos::ParameterList& fsitimeparams = manipulate_fsi_time_params(fs3idyn);

  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    {
      fsi_ = Teuchos::rcp(new FSI::MonolithicFluidSplit(comm_, fsitimeparams));
      break;
    }
    case fsi_iter_monolithicstructuresplit:
    {
      fsi_ = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm_, fsitimeparams));
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown FSI coupling algorithm");
      break;
    }
  }

  //---------------------------------------------------------------------
  // create instances for fluid- and structure-based scalar transport
  // solver and arrange them in combined vector
  //---------------------------------------------------------------------
  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3idyn.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
  const int linsolver2number = fs3idyn.get<int>("LINEAR_SOLVER2");

  // check if the linear solver has a valid solver number
  if (linsolver1number == (-1))
    FOUR_C_THROW(
        "no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER2 in "
        "FS3I DYNAMIC to a valid number!");
  if (linsolver2number == (-1))
    FOUR_C_THROW(
        "no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 "
        "in FS3I DYNAMIC to a valid number!");

  fluidscatra_ = Teuchos::rcp(
      new Adapter::ScaTraBaseAlgorithm(fs3idyn, problem->scalar_transport_dynamic_params(),
          problem->SolverParams(linsolver1number), "scatra1", true));
  fluidscatra_->init();
  fluidscatra_->ScaTraField()->set_number_of_dof_set_displacement(1);
  fluidscatra_->ScaTraField()->set_number_of_dof_set_velocity(1);
  fluidscatra_->ScaTraField()->set_number_of_dof_set_wall_shear_stress(1);

  structscatra_ = Teuchos::rcp(
      new Adapter::ScaTraBaseAlgorithm(fs3idyn, problem->scalar_transport_dynamic_params(),
          problem->SolverParams(linsolver2number), "scatra2", true));
  structscatra_->init();
  structscatra_->ScaTraField()->set_number_of_dof_set_displacement(1);
  structscatra_->ScaTraField()->set_number_of_dof_set_velocity(1);
  structscatra_->ScaTraField()->set_number_of_dof_set_wall_shear_stress(1);

  scatravec_.push_back(fluidscatra_);
  scatravec_.push_back(structscatra_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::setup()
{
  // call setup in base class
  FS3I::FS3IBase::setup();

  // setup structure scatra
  structscatra_->setup();

  // setup fluid scatra
  fluidscatra_->setup();

  //---------------------------------------------------------------------
  // check existence of scatra coupling conditions for both
  // discretizations and definition of the permeability coefficient
  //---------------------------------------------------------------------
  CheckFS3IInputs();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Adapter::MortarVolCoupl> FS3I::PartFS3I::create_vol_mortar_object(
    Teuchos::RCP<Core::FE::Discretization> masterdis,
    Teuchos::RCP<Core::FE::Discretization> slavedis)
{
  // copy conditions
  // this is actually only needed for copying TRANSPORT DIRICHLET/NEUMANN CONDITIONS
  // as standard DIRICHLET/NEUMANN CONDITIONS
  ScaTra::ScatraFluidCloneStrategy clonestrategy;
  const auto conditions_to_copy = clonestrategy.conditions_to_copy();
  Core::FE::DiscretizationCreatorBase creator;
  creator.CopyConditions(*slavedis, *slavedis, conditions_to_copy);

  // first call fill_complete for single discretizations.
  // This way the physical dofs are numbered successively
  masterdis->fill_complete();
  slavedis->fill_complete();

  // build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_scatra = slavedis->NumDof(0, slavedis->lRowNode(0));
  const int ndofperelement_scatra = 0;
  const int ndofpernode_struct = masterdis->NumDof(0, masterdis->lRowNode(0));
  const int ndofperelement_struct = 0;

  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
  dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
      ndofpernode_scatra, ndofperelement_scatra, 0, true));
  if (masterdis->AddDofSet(dofsetaux) != 1) FOUR_C_THROW("unexpected dof sets in structure field");
  dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
      ndofpernode_struct, ndofperelement_struct, 0, true));
  if (slavedis->AddDofSet(dofsetaux) != 1) FOUR_C_THROW("unexpected dof sets in scatra field");

  // call assign_degrees_of_freedom also for auxiliary dofsets
  // note: the order of fill_complete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. scatra dofs
  // 3. structure auxiliary dofs
  // 4. scatra auxiliary dofs
  masterdis->fill_complete(true, false, false);
  slavedis->fill_complete(true, false, false);


  // Scheme: non matching meshes --> volumetric mortar coupling...
  Teuchos::RCP<Core::Adapter::MortarVolCoupl> volume_coupling_object =
      Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

  // setup projection matrices (use default material strategy)
  volume_coupling_object->init(Global::Problem::Instance()->NDim(), masterdis, slavedis);
  Teuchos::ParameterList binning_params = Global::Problem::Instance()->binning_strategy_params();
  Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
      "spatial_approximation_type", Global::Problem::Instance()->spatial_approximation_type(),
      binning_params);
  volume_coupling_object->Redistribute(
      binning_params, Global::Problem::Instance()->OutputControlFile());
  volume_coupling_object->setup(Global::Problem::Instance()->VolmortarParams());

  return volume_coupling_object;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList& FS3I::PartFS3I::manipulate_fsi_time_params(
    const Teuchos::ParameterList& fs3idyn) const
{
  // NOTE: we can not do this in the AC-fs3i class were it would belong,
  // since overloading a function inside the constructor does not work :(

  Teuchos::ParameterList& timeparams = *(new Teuchos::ParameterList(fs3idyn));

  const int fsisubcycles = fs3idyn.sublist("AC").get<int>("FSI_STEPS_PER_SCATRA_STEP");

  if (fsisubcycles != 1)  // if we have subcycling for ac_fsi
  {
    timeparams.set<double>("TIMESTEP", fs3idyn.get<double>("TIMESTEP") / (double)fsisubcycles);
  }
  return timeparams;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::read_restart()
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  const int restart = Global::Problem::Instance()->restart();

  if (restart)
  {
    const Teuchos::ParameterList& fs3idynac = Global::Problem::Instance()->FS3IDynamicParams();
    const bool restartfrompartfsi =
        Core::UTILS::IntegralValue<int>(fs3idynac, "RESTART_FROM_PART_FSI");

    if (not restartfrompartfsi)  // standard restart
    {
      fsi_->read_restart(restart);

      for (unsigned i = 0; i < scatravec_.size(); ++i)
      {
        Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
        currscatra->ScaTraField()->read_restart(restart);
      }
    }
    else  // we do not want to read the scatras values and the lagrange multiplyer, since we start
          // from a partitioned FSI
    {
      fsi_->read_restart(restart);

      // we need to set time and step in scatra to have it matching with the fsi ones
      for (unsigned i = 0; i < scatravec_.size(); ++i)
      {
        Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
        currscatra->ScaTraField()->SetTimeStep(
            fsi_->fluid_field()->Time(), fsi_->fluid_field()->Step());
      }
    }

    time_ = fsi_->fluid_field()->Time();
    step_ = fsi_->fluid_field()->Step();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetupSystem()
{
  // now do the coupling setup and create the combined dofmap
  fsi_->SetupSystem();

  /*----------------------------------------------------------------------*/
  /*                  General set up for scalar fields                    */
  /*----------------------------------------------------------------------*/

  // create map extractors needed for scatra condition coupling
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
    Teuchos::RCP<Core::FE::Discretization> currdis = currscatra->ScaTraField()->discretization();
    const int numscal = currscatra->ScaTraField()->NumScal();
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> mapex =
        Teuchos::rcp(new Core::LinAlg::MultiMapExtractor());
    Core::Conditions::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(
        new Core::Conditions::NDimConditionSelector(*currdis, "ScaTraCoupling", 0, numscal)));
    mcs.SetupExtractor(*currdis, *currdis->dof_row_map(), *mapex);
    scatrafieldexvec_.push_back(mapex);
  }

  scatracoup_->setup_condition_coupling(*(scatravec_[0]->ScaTraField()->discretization()),
      scatrafieldexvec_[0]->Map(1), *(scatravec_[1]->ScaTraField()->discretization()),
      scatrafieldexvec_[1]->Map(1), "ScaTraCoupling",
      scatravec_[0]
          ->ScaTraField()
          ->NumScal());  // we assume here that both discretisation have the same number of scalars

  // create map extractor for coupled scatra fields
  // the second field (currently structure) is always split
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;

  // In the limiting case of an infinite permeability of the interface between
  // different scatra fields, the concentrations on both sides of the interface are
  // constrained to be equal. In this case, we keep the fluid scatra dofs at the
  // interface as unknowns in the overall system, whereas the structure scatra
  // dofs are condensed (cf. "structuresplit" in a monolithic FSI
  // system). Otherwise, both concentrations are kept in the overall system
  // and the equality of fluxes is considered explicitly.
  if (infperm_)
  {
    maps.push_back(scatrafieldexvec_[0]->FullMap());
    maps.push_back(scatrafieldexvec_[1]->Map(0));
  }
  else
  {
    maps.push_back(scatrafieldexvec_[0]->FullMap());
    maps.push_back(scatrafieldexvec_[1]->FullMap());
  }
  Teuchos::RCP<Epetra_Map> fullmap = Core::LinAlg::MultiMapExtractor::MergeMaps(maps);
  scatraglobalex_->setup(*fullmap, maps);

  // create coupling vectors and matrices (only needed for finite surface permeabilities)
  if (!infperm_)
  {
    for (unsigned i = 0; i < scatravec_.size(); ++i)
    {
      Teuchos::RCP<Epetra_Vector> scatracoupforce =
          Teuchos::rcp(new Epetra_Vector(*(scatraglobalex_->Map(i)), true));
      scatracoupforce_.push_back(scatracoupforce);

      Teuchos::RCP<Core::LinAlg::SparseMatrix> scatracoupmat =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(scatraglobalex_->Map(i)), 27, false, true));
      scatracoupmat_.push_back(scatracoupmat);

      const Epetra_Map* dofrowmap = scatravec_[i]->ScaTraField()->discretization()->dof_row_map();
      Teuchos::RCP<Epetra_Vector> zeros = Core::LinAlg::CreateVector(*dofrowmap, true);
      scatrazeros_.push_back(zeros);
    }
  }

  // create scatra block matrix
  scatrasystemmatrix_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *scatraglobalex_, *scatraglobalex_, 27, false, true));

  // create scatra rhs vector
  scatrarhs_ = Teuchos::rcp(new Epetra_Vector(*scatraglobalex_->FullMap(), true));

  // create scatra increment vector
  scatraincrement_ = Teuchos::rcp(new Epetra_Vector(*scatraglobalex_->FullMap(), true));

  // check whether potential Dirichlet conditions at scatra interface are
  // defined for both discretizations
  check_interface_dirichlet_bc();

  // scatra solver
  Teuchos::RCP<Core::FE::Discretization> firstscatradis =
      (scatravec_[0])->ScaTraField()->discretization();
#ifdef SCATRABLOCKMATRIXMERGE
  Teuchos::RCP<Teuchos::ParameterList> scatrasolvparams = Teuchos::rcp(new Teuchos::ParameterList);
  Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
      "SOLVER", Core::LinearSolver::SolverType::umfpack, scatrasolvparams);
  scatrasolver_ = Teuchos::rcp(new Core::LinAlg::Solver(scatrasolvparams, firstscatradis->Comm()));
#else
  const Teuchos::ParameterList& fs3idyn = Global::Problem::Instance()->FS3IDynamicParams();
  // get solver number used for fs3i
  const int linsolvernumber = fs3idyn.get<int>("COUPLED_LINEAR_SOLVER");
  // check if LOMA solvers has a valid number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for FS3I problems. Please set COUPLED_LINEAR_SOLVER in FS3I "
        "DYNAMIC to a valid number!");

  const Teuchos::ParameterList& coupledscatrasolvparams =
      Global::Problem::Instance()->SolverParams(linsolvernumber);

  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(coupledscatrasolvparams, "SOLVER");

  if (solvertype != Core::LinearSolver::SolverType::belos)
    FOUR_C_THROW("Iterative solver expected");

  const auto azprectype = Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(
      coupledscatrasolvparams, "AZPREC");
  if (azprectype != Core::LinearSolver::PreconditionerType::block_gauss_seidel_2x2)
    FOUR_C_THROW("Block Gauss-Seidel preconditioner expected");

  // use coupled scatra solver object
  scatrasolver_ = Teuchos::rcp(new Core::LinAlg::Solver(coupledscatrasolvparams,
      firstscatradis->Comm(), Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY")));

  // get the solver number used for fluid ScalarTransport solver
  const int linsolver1number = fs3idyn.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
  const int linsolver2number = fs3idyn.get<int>("LINEAR_SOLVER2");

  // check if the linear solver has a valid solver number
  if (linsolver1number == (-1))
    FOUR_C_THROW(
        "no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER1 in "
        "FS3I DYNAMIC to a valid number!");
  if (linsolver2number == (-1))
    FOUR_C_THROW(
        "no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 "
        "in FS3I DYNAMIC to a valid number!");

  scatrasolver_->put_solver_params_to_sub_params("Inverse1",
      Global::Problem::Instance()->SolverParams(linsolver1number),
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));
  scatrasolver_->put_solver_params_to_sub_params("Inverse2",
      Global::Problem::Instance()->SolverParams(linsolver2number),
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));

  (scatravec_[0])
      ->ScaTraField()
      ->discretization()
      ->compute_null_space_if_necessary(scatrasolver_->Params().sublist("Inverse1"));
  (scatravec_[1])
      ->ScaTraField()
      ->discretization()
      ->compute_null_space_if_necessary(scatrasolver_->Params().sublist("Inverse2"));
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::TestResults(const Epetra_Comm& comm)
{
  Global::Problem::Instance()->AddFieldTest(fsi_->fluid_field()->CreateFieldTest());
  Global::Problem::Instance()->AddFieldTest(fsi_->ale_field()->CreateFieldTest());
  Global::Problem::Instance()->AddFieldTest(fsi_->structure_field()->CreateFieldTest());

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    Global::Problem::Instance()->AddFieldTest(scatra->create_sca_tra_field_test());
  }
  Global::Problem::Instance()->TestAll(comm);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetFSISolution()
{
  // we clear every state, including the states of the secondary dof sets
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    scatravec_[i]->ScaTraField()->discretization()->ClearState(true);
    // we have to manually clear this since this can not be saved directly in the
    // primary dof set (because it is cleared in between)
    scatravec_[i]->ScaTraField()->clear_external_concentrations();
  }

  set_mesh_disp();
  set_velocity_fields();
  set_wall_shear_stresses();
  set_membrane_concentration();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::set_struct_scatra_solution() const
{
  fsi_->structure_field()->discretization()->set_state(
      1, "scalarfield", structure_scalar_to_structure(scatravec_[1]->ScaTraField()->Phinp()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::set_mesh_disp() const
{
  // fluid field
  Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> fluidscatra = scatravec_[0];
  fluidscatra->ScaTraField()->ApplyMeshMovement(FluidToFluidScalar(fsi_->fluid_field()->Dispnp()));

  // structure field
  Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> structscatra = scatravec_[1];
  structscatra->ScaTraField()->ApplyMeshMovement(
      structure_to_structure_scalar(fsi_->structure_field()->Dispnp()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::set_velocity_fields() const
{
  std::vector<Teuchos::RCP<const Epetra_Vector>> convel;
  std::vector<Teuchos::RCP<const Epetra_Vector>> vel;
  ExtractVel(convel, vel);

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->set_velocity_field(vol_mortar_master_to_slavei(i, convel[i]),
        Teuchos::null, vol_mortar_master_to_slavei(i, vel[i]), Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector>>& convel,
    std::vector<Teuchos::RCP<const Epetra_Vector>>& vel) const
{
  // extract fluid velocities

  switch (fsi_->fluid_field()->TimIntScheme())
  {
    case Inpar::FLUID::timeint_afgenalpha:
    {
      Teuchos::RCP<Epetra_Vector> fluidconvel =
          Teuchos::rcp(new Epetra_Vector(*(fsi_->fluid_field()->Velaf())));
      vel.push_back(fluidconvel);
      // now subtract the grid velocity
      fluidconvel->Update(-1.0, *(fsi_->fluid_field()->GridVel()), 1.0);
      convel.push_back(fluidconvel);
    }
    break;
    case Inpar::FLUID::timeint_one_step_theta:
    {
      convel.push_back(fsi_->fluid_field()->ConvectiveVel());
      vel.push_back(fsi_->fluid_field()->Velnp());
    }
    break;
    default:
      FOUR_C_THROW("Time integration scheme not supported");
      break;
  }

  // extract structure velocities and accelerations

  Teuchos::RCP<Epetra_Vector> velocity =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->Velnp())));
  vel.push_back(velocity);
  // structure ScaTra: velocity and grid velocity are identical!
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(velocity->Map(), true));
  convel.push_back(zeros);
}

/*----------------------------------------------------------------------*
 |  Set wall shear stresses                                  Thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFS3I::set_wall_shear_stresses() const
{
  std::vector<Teuchos::RCP<const Epetra_Vector>> wss;
  ExtractWSS(wss);

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->set_wall_shear_stresses(vol_mortar_master_to_slavei(i, wss[i]));
  }
}

/*----------------------------------------------------------------------*
 |  Extract wall shear stresses                              Thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFS3I::ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector>>& wss) const
{
  // ############ Fluid Field ###############

  Teuchos::RCP<Adapter::FluidFSI> fluid =
      Teuchos::rcp_dynamic_cast<Adapter::FluidFSI>(fsi_->fluid_field());
  if (fluid == Teuchos::null) FOUR_C_THROW("Dynamic cast to Adapter::FluidFSI failed!");

  Teuchos::RCP<Epetra_Vector> WallShearStress = fluid->calculate_wall_shear_stresses();

  if (Core::UTILS::IntegralValue<Inpar::FLUID::WSSType>(
          Global::Problem::Instance()->FluidDynamicParams(), "WSS_TYPE") !=
      Inpar::FLUID::wss_standard)
    FOUR_C_THROW("WSS_TYPE not supported for FS3I!");

  wss.push_back(WallShearStress);

  // ############ Structure Field ###############

  // extract FSI-Interface from fluid field
  WallShearStress = fsi_->fluid_field()->Interface()->ExtractFSICondVector(WallShearStress);

  // replace global fluid interface dofs through structure interface dofs
  WallShearStress = fsi_->fluid_to_struct(WallShearStress);

  // insert structure interface entries into vector with full structure length
  Teuchos::RCP<Epetra_Vector> structure =
      Core::LinAlg::CreateVector(*(fsi_->structure_field()->Interface()->FullMap()), true);

  // Parameter int block of function InsertVector: (0: inner dofs of structure, 1: interface dofs of
  // structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
  fsi_->structure_field()->Interface()->InsertVector(WallShearStress, 1, structure);
  wss.push_back(structure);
}

/*----------------------------------------------------------------------*
 |  transport quantity from fluid to fluid-scalar            Thon 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::FluidToFluidScalar(
    const Teuchos::RCP<const Epetra_Vector> fluidvector) const
{
  return vol_mortar_master_to_slavei(0, fluidvector);
}

/*----------------------------------------------------------------------*
 |  transport quantity from fluid-scalar to fluid            Thon 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::FluidScalarToFluid(
    const Teuchos::RCP<const Epetra_Vector> fluidscalarvector) const
{
  return vol_mortar_slave_to_masteri(0, fluidscalarvector);
}

/*----------------------------------------------------------------------*
 |  transport quantity from structure to structure-scalar    Thon 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::structure_to_structure_scalar(
    const Teuchos::RCP<const Epetra_Vector> structurevector) const
{
  return vol_mortar_master_to_slavei(1, structurevector);
}

/*----------------------------------------------------------------------*
 |  transport quantity from structure-scalar to structure    Thon 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::structure_scalar_to_structure(
    const Teuchos::RCP<const Epetra_Vector> structurescalavector) const
{
  return vol_mortar_slave_to_masteri(1, structurescalavector);
}

/*-------------------------------------------------------------------------------------*
 |  transport quantity from i-th volmortar master to i-th volmortar slave   Thon 08/16 |
 *-------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::vol_mortar_master_to_slavei(
    const int i, const Teuchos::RCP<const Epetra_Vector> mastervector) const
{
  switch (volume_fieldcouplings_[i])
  {
    case Inpar::FS3I::coupling_match:
      return mastervector;
      break;
    case Inpar::FS3I::coupling_nonmatch:
      return volume_coupling_objects_[i]->apply_vector_mapping21(mastervector);
      break;
    default:
      FOUR_C_THROW("unknown field coupling type");
      return Teuchos::null;
      break;
  }
}

/*-------------------------------------------------------------------------------------*
 |  transport quantity from i-th volmortar slave to i-th volmortar master   Thon 08/16 |
 *-------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::vol_mortar_slave_to_masteri(
    const int i, const Teuchos::RCP<const Epetra_Vector> slavevector) const
{
  switch (volume_fieldcouplings_[i])
  {
    case Inpar::FS3I::coupling_match:
      return slavevector;
      break;
    case Inpar::FS3I::coupling_nonmatch:
      return volume_coupling_objects_[i]->apply_vector_mapping12(slavevector);
      break;
    default:
      FOUR_C_THROW("unknown field coupling type");
      return Teuchos::null;
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
