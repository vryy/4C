/*-----------------------------------------------------------*/
/*!
\file ad_str_structure_new.cpp

\brief Adapter for the new structural time integration framework.

\maintainer Michael Hiermeier

\date Sep 2, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "ad_str_structure_new.H"
#include "ad_str_timint_adaptive.H"
#include "ad_str_fsi_timint_adaptive.H"
#include "ad_str_constr_merged.H"
#include "ad_str_wrapper.H"
#include "ad_str_lung.H"
#include "ad_str_redairway.H"
#include "ad_str_fsi_crack.H"
#include "ad_str_fpsiwrapper.H"
#include "ad_str_fsiwrapper_immersed.H"
#include "ad_str_invana.H"

#include "../drt_structure/strtimada_create.H"

#include "../drt_crack/InsertCohesiveElements.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_structure_new/str_timint_factory.H"
#include "../drt_structure_new/str_solver_factory.H"
#include "../drt_structure_new/str_timint_base.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_mat/matpar_bundle.H"

#include "../drt_inpar/inpar_crack.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/drt_validparameters.H"

#include "../drt_so3/so_sh8p8.H"
#include "../drt_so3/so3_ssn_plast_eletypes.H"
#include "../drt_so3/so3_ssn_plast_sosh8.H"
#include "../drt_so3/so3_ssn_plast_sosh18.H"
#include "../drt_so3/so_hex8fbar.H"
#include "../drt_s8/shell8.H"

#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3k.H"

#include "../solver_nonlin_nox/nox_nln_group.H"
#include "../solver_nonlin_nox/nox_nln_group_prepostoperator.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ADAPTER::StructureBaseAlgorithmNew::StructureBaseAlgorithmNew()
    : str_wrapper_(Teuchos::null),
      prbdyn_(Teuchos::null),
      sdyn_(Teuchos::null),
      isinit_(false),
      issetup_(false)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::Init(
    const Teuchos::ParameterList& prbdyn,
    Teuchos::ParameterList& sdyn,
    Teuchos::RCP<DRT::Discretization> actdis)
{
  issetup_ = false;

  // save a pointer of the input variables as class variables
  prbdyn_ = Teuchos::rcpFromRef<const Teuchos::ParameterList>(prbdyn);
  sdyn_ = Teuchos::rcpFromRef<Teuchos::ParameterList>(sdyn);
  actdis_ = actdis;

  isinit_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::Setup()
{
  if (not IsInit())
    dserror("You have to call Init() first!");

  // major switch to different time integrators
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(*sdyn_,"DYNAMICTYP"))
  {
  case INPAR::STR::dyna_statics :
  case INPAR::STR::dyna_genalpha :
  case INPAR::STR::dyna_genalpha_liegroup :
  case INPAR::STR::dyna_onesteptheta :
  case INPAR::STR::dyna_gemm :
  case INPAR::STR::dyna_expleuler:
  case INPAR::STR::dyna_centrdiff :
  case INPAR::STR::dyna_ab2 :
  case INPAR::STR::dyna_euma :
  case INPAR::STR::dyna_euimsto :
    SetupTimInt();  // <-- here is the show
    break;
  default :
    dserror("Unknown time integration scheme '%s'", sdyn_->
        get<std::string>("DYNAMICTYP").c_str());
    break;
  }

  issetup_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::SetupTimInt()
{
  if (not IsInit())
    dserror("You have to call Init() first!");

  // get the problem instance
  DRT::Problem* problem = DRT::Problem::Instance();
  // what's the current problem type?
  PROBLEM_TYP probtype = problem->ProblemType();
  // get the restart step
  const int restart = problem->Restart();

  // ---------------------------------------------------------------------------
  // Define, initialize and start the timer
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::Time> t
    = Teuchos::TimeMonitor::getNewTimer(
        "ADAPTER::StructureTimIntBaseAlgorithm::SetupStructure");
  Teuchos::TimeMonitor monitor(*t);

  // ---------------------------------------------------------------------------
  // Here we read the discretization at the current
  // time step from restart files
  // ---------------------------------------------------------------------------
  if ( restart and probtype == prb_crack )
  {
    IO::DiscretizationReader reader(actdis_, restart);
    reader.ReadMesh(restart);
  }
  // set degrees of freedom in the discretization
  else if (not actdis_->Filled() || not actdis_->HaveDofs())
    actdis_->FillComplete();

  // ---------------------------------------------------------------------------
  // Setup a model type set by checking
  // the different conditions
  // ---------------------------------------------------------------------------
  // define and initial with default value
  Teuchos::RCP<std::set<enum INPAR::STR::ModelType> > modeltypes =
      Teuchos::rcp(new std::set<enum INPAR::STR::ModelType>());
  modeltypes->insert(INPAR::STR::model_structure);
  SetModelTypes(*modeltypes);

  // ---------------------------------------------------------------------------
  // Setup a element technology set by checking
  // the elements of the discretization
  // ---------------------------------------------------------------------------
  Teuchos::RCP<std::set<enum INPAR::STR::EleTech> > eletechs =
      Teuchos::rcp(new std::set<enum INPAR::STR::EleTech>());
  DetectElementTechnologies(*eletechs);

  // ---------------------------------------------------------------------------
  // Setup the parameter lists for structural
  // time integration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> ioflags
    = Teuchos::rcp(new Teuchos::ParameterList(problem->IOParams()));
  Teuchos::RCP<Teuchos::ParameterList> taflags
    = Teuchos::rcp(new Teuchos::ParameterList(sdyn_->sublist("TIMEADAPTIVITY")));
  Teuchos::RCP<Teuchos::ParameterList> xparams
    = Teuchos::rcp(new Teuchos::ParameterList());
  SetParams(*ioflags,*xparams,*taflags);

  // ---------------------------------------------------------------------------
  // Setup and create model specific linear solvers
  // ---------------------------------------------------------------------------
  Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver> > > linsolvers =
      STR::SOLVER::BuildLinSolvers(*modeltypes,*sdyn_,*actdis_);

  // ---------------------------------------------------------------------------
  // Checks in case of multi-scale simulations
  // ---------------------------------------------------------------------------
  {
    // make sure we IMR-like generalised-alpha requested for multi-scale
    // simulations
    Teuchos::RCP<MAT::PAR::Bundle> materials = problem->Materials();
    for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials->Map()->begin();
         i!=materials->Map()->end();
         ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> mat = i->second;
      if (mat->Type() == INPAR::MAT::m_struct_multiscale)
      {
        if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(*sdyn_, "DYNAMICTYP") != INPAR::STR::dyna_genalpha)
          dserror("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha");
        else if (DRT::INPUT::IntegralValue<INPAR::STR::MidAverageEnum>(sdyn_->sublist("GENALPHA"), "GENAVG") != INPAR::STR::midavg_trlike)
          dserror("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha with GENAVG=TrLike");
        break;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Add cohesive elements in case of crack
  // propagation simulations
  // ---------------------------------------------------------------------------
  {
    const Teuchos::ParameterList& crackparam = DRT::Problem::Instance()->CrackParams();
    if (DRT::INPUT::IntegralValue<INPAR::CRACK::crackModel>(crackparam,"CRACK_MODEL")
          == INPAR::CRACK::crack_cohesive)
    {
      //if( (not DRT::Problem::Instance()->Restart()) and DRT::Problem::Instance()->ProblemType() == prb_structure )
      {
        Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
        DRT::CRACK::InsertCohesiveElements isp( structdis );
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Create context for output and restart
  // ---------------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis_->Writer();
  if (DRT::INPUT::IntegralValue<int>(*ioflags,"OUTPUT_BIN"))
  {
    output->WriteMesh(0, 0.0);
  }

  // ---------------------------------------------------------------------------
  // initialize/setup the input/output data container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio =
      Teuchos::rcp(new STR::TIMINT::BaseDataIO());
  dataio->Init(*ioflags,*sdyn_,*xparams,output);
  dataio->Setup();

  // ---------------------------------------------------------------------------
  // initialize/setup the structural dynamics data
  // container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn =
      STR::TIMINT::BuildDataSDyn(*sdyn_);
  datasdyn->Init(actdis_,*sdyn_,*xparams,modeltypes,eletechs,linsolvers);
  datasdyn->Setup();

  // ---------------------------------------------------------------------------
  // initialize/setup the global state data container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate =
      Teuchos::rcp(new STR::TIMINT::BaseDataGlobalState());
  dataglobalstate->Init(actdis_,*sdyn_,datasdyn);
  dataglobalstate->Setup();

  // ---------------------------------------------------------------------------
  // in case of non-additive rotation (pseudo-)vector DOFs:
  // ---------------------------------------------------------------------------
  if (eletechs->find(INPAR::STR::eletech_rotvec) != eletechs->end())
  {
    // -------------------------------------------------------------------------
    // setup the map extractor for split additive<->rotvec DOFs
    // -------------------------------------------------------------------------
    dataglobalstate->SetupRotVecMapExtractor();

    // -------------------------------------------------------------------------
    // set the RotVecUpdater as new precomputeX operator for the nox nln group
    // -------------------------------------------------------------------------
    Teuchos::ParameterList& p_grp_opt =
        datasdyn->GetMutableNoxParams().sublist("Group Options");
    // Get the current map. If there is no map, return a new empty one. (reference)
    NOX::NLN::GROUP::PrePostOperator::Map& prepostgroup_map =
        NOX::NLN::GROUP::PrePostOp::GetMutableMap(p_grp_opt);
    // create the new rotation vector update pre/post operator
    Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> prepostrotvec_ptr =
        Teuchos::rcp(new NOX::NLN::GROUP::PrePostOp::TIMINT::RotVecUpdater(
            dataglobalstate));
    // insert/replace the old pointer in the map
    prepostgroup_map[NOX::NLN::GROUP::prepost_rotvecupdate] = prepostrotvec_ptr;
  }

  // ---------------------------------------------------------------------------
  // Build time integrator
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TIMINT::Base> ti_strategy =
      STR::TIMINT::BuildStrategy(*sdyn_);
  ti_strategy->Init(dataio,datasdyn,dataglobalstate);
  /* In the restart case, we Setup the structural time integration after the
   * discretization has been redistributed. See STR::TIMINT::Base::ReadRestart()
   * for more information.                                     hiermeier 05/16*/
  if (not restart)
    ti_strategy->Setup();


  // ---------------------------------------------------------------------------
  // Create wrapper for the time integration strategy
  // ---------------------------------------------------------------------------
  SetStructureWrapper(*ioflags,*sdyn_,*xparams,*taflags,ti_strategy);

  // see you
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::SetModelTypes(
    std::set<enum INPAR::STR::ModelType>& modeltypes) const
{
  if (not IsInit())
    dserror("You have to call Init() first!");

  // ---------------------------------------------------------------------------
  // check for meshtying and contact conditions
  // ---------------------------------------------------------------------------
  // --- contact conditions
  std::vector<DRT::Condition*> ccond(0);
  actdis_->GetCondition("Contact",ccond);
  if (ccond.size())
    modeltypes.insert(INPAR::STR::model_contact);
  // --- meshtying conditions
  std::vector<DRT::Condition*> mtcond(0);
  actdis_->GetCondition("Mortar", mtcond);
  if (mtcond.size())
    modeltypes.insert(INPAR::STR::model_meshtying);

  // check for 0D cardiovascular conditions
  // ---------------------------------------------------------------------------
  std::vector<DRT::Condition*> cardiovasc0dcond_windkesselonly(0);
  std::vector<DRT::Condition*> cardiovasc0dcond_arterialproxdist(0);
  std::vector<DRT::Condition*> cardiovasc0dcond_arterialvenoussyspulcoupled(0);
  actdis_->GetCondition("Cardiovascular0DWindkesselOnlyStructureCond",cardiovasc0dcond_windkesselonly);
  actdis_->GetCondition("Cardiovascular0DArterialProxDistStructureCond",
      cardiovasc0dcond_arterialproxdist);
  actdis_->GetCondition("Cardiovascular0DArterialVenousSysPulCoupledStructureCond",
      cardiovasc0dcond_arterialvenoussyspulcoupled);
  if (cardiovasc0dcond_windkesselonly.size() or
      cardiovasc0dcond_arterialproxdist.size() or
      cardiovasc0dcond_arterialvenoussyspulcoupled.size())
    modeltypes.insert(INPAR::STR::model_cardiovascular0d);

  // ---------------------------------------------------------------------------
  // check for constraint conditions
  // ---------------------------------------------------------------------------
  bool have_lag_constraint = false;
  bool have_pen_constraint = false;
  // --- enforcement by Lagrange multiplier
  std::vector<DRT::Condition*> lagcond_volconstr3d(0);
  std::vector<DRT::Condition*> lagcond_areaconstr3d(0);
  std::vector<DRT::Condition*> lagcond_areaconstr2d(0);
  std::vector<DRT::Condition*> lagcond_mpconline2d(0);
  std::vector<DRT::Condition*> lagcond_mpconplane3d(0);
  std::vector<DRT::Condition*> lagcond_mpcnormcomp3d(0);
  actdis_->GetCondition("VolumeConstraint_3D",lagcond_volconstr3d);
  actdis_->GetCondition("AreaConstraint_3D",lagcond_areaconstr3d);
  actdis_->GetCondition("AreaConstraint_2D",lagcond_areaconstr2d);
  actdis_->GetCondition("MPC_NodeOnLine_2D",lagcond_mpconline2d);
  actdis_->GetCondition("MPC_NodeOnPlane_3D",lagcond_mpconplane3d);
  actdis_->GetCondition("MPC_NormalComponent_3D",lagcond_mpcnormcomp3d);
  if (
         lagcond_volconstr3d.size()  or
         lagcond_areaconstr3d.size() or
         lagcond_areaconstr2d.size() or
         lagcond_mpconline2d.size()  or
         lagcond_mpconplane3d.size() or
         lagcond_mpcnormcomp3d.size()
      )
    have_lag_constraint = true;
  // --- enforcement by penalty law
  std::vector<DRT::Condition*> pencond_volconstr3d(0);
  std::vector<DRT::Condition*> pencond_areaconstr3d(0);
  std::vector<DRT::Condition*> pencond_mpcnormcomp3d(0);
  actdis_->GetCondition("VolumeConstraint_3D_Pen",pencond_volconstr3d);
  actdis_->GetCondition("AreaConstraint_3D_Pen",pencond_areaconstr3d);
  actdis_->GetCondition("MPC_NormalComponent_3D_Pen",pencond_mpcnormcomp3d);
  if (
         pencond_volconstr3d.size() or
         pencond_areaconstr3d.size() or
         pencond_mpcnormcomp3d.size()
      )
    have_pen_constraint = true;
  if (have_lag_constraint or have_pen_constraint)
    modeltypes.insert(INPAR::STR::model_lag_pen_constraint);

  // ---------------------------------------------------------------------------
  // check for spring dashpot conditions
  // ---------------------------------------------------------------------------
  std::vector<DRT::Condition*> sdp_cond(0);
  actdis_->GetCondition("RobinSpringDashpot",sdp_cond);
  if (sdp_cond.size())
    modeltypes.insert(INPAR::STR::model_springdashpot);
  // ---------------------------------------------------------------------------
  // check for coupled problems
  // ---------------------------------------------------------------------------
  // get the problem instance
  DRT::Problem* problem = DRT::Problem::Instance();
  // what's the current problem type?
  PROBLEM_TYP probtype = problem->ProblemType();
  switch (probtype)
  {
    case prb_fsi:
    case prb_immersed_fsi:
    case prb_immersed_ale_fsi:
    case prb_fsi_redmodels:
    case prb_fsi_lung:
    case prb_gas_fsi:
    case prb_ac_fsi:
    case prb_biofilm_fsi:
    case prb_thermo_fsi:
    case prb_fsi_xfem:
    case prb_fsi_crack:
    {
      if (prbdyn_->INVALID_TEMPLATE_QUALIFIER
          isType<Teuchos::RCP<STR::MODELEVALUATOR::Generic> > ("Partitioned Coupling Model"))
      {
        const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& coupling_model_ptr =
          prbdyn_->INVALID_TEMPLATE_QUALIFIER
          get<Teuchos::RCP<STR::MODELEVALUATOR::Generic> >("Partitioned Coupling Model");
        if (coupling_model_ptr.is_null())
          dserror("The partitioned coupling model pointer is not allowed to be Teuchos::null!");
        // set the model type
        modeltypes.insert(INPAR::STR::model_partitioned_coupling);
        // copy the coupling model object pointer into the (temporal) sdyn parameter list
        sdyn_->set<Teuchos::RCP<STR::MODELEVALUATOR::Generic> >("Partitioned Coupling Model",
            coupling_model_ptr);
      }
      break;
    }
    default:
      // do nothing
      break;
  } // switch (probtype)

  // ---------------------------------------------------------------------------
  // check for beam interactions (either contact or potential-based)
  // ---------------------------------------------------------------------------
  // get beam contact strategy since there are no conditions for beam contact
  const Teuchos::ParameterList& beamcontact =
     DRT::Problem::Instance()->BeamContactParams();
  INPAR::BEAMCONTACT::Strategy strategy =
      DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(beamcontact,
          "BEAMS_STRATEGY");

  // conditions for potential-based beam interaction
  std::vector<DRT::Condition*> beampotconditions(0);
  actdis_->GetCondition("BeamPotentialLineCharge",beampotconditions);

  if (strategy != INPAR::BEAMCONTACT::bstr_none or beampotconditions.size())
    modeltypes.insert(INPAR::STR::model_beam_interaction);

  // ---------------------------------------------------------------------------
  // check for brownian dynamics
  // ---------------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(), "STATMECHPROB"))
    modeltypes.insert(INPAR::STR::model_browniandyn);

  // ---------------------------------------------------------------------------
  // check for crosslinking in biopolymer networks
  // ---------------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(), "CROSSLINKER"))
    modeltypes.insert(INPAR::STR::model_crosslinking);

  // hopefully we haven't forgotten anything
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::DetectElementTechnologies(
    std::set<enum INPAR::STR::EleTech>& eletechs) const
{
  int isplasticity_local = 0;
  int isplasticity_global = 0;

  int iseas_local = 0;
  int iseas_global = 0;

  int isfbar_local = 0;
  int isfbar_global = 0;

  int ispressure_local = 0;
  int ispressure_global = 0;

  int isrotvec_local = 0;
  int isrotvec_global = 0;

  for (int i=0;i<actdis_->NumMyRowElements();++i)
  {
    DRT::Element* actele = actdis_->lRowElement(i);
    // Detect plasticity -------------------------------------------------------
    if (actele->ElementType() == DRT::ELEMENTS::So_hex8PlastType::Instance() or
        actele->ElementType() == DRT::ELEMENTS::So_hex27PlastType::Instance() or
        actele->ElementType() == DRT::ELEMENTS::So_sh8PlastType::Instance() or
        actele->ElementType() == DRT::ELEMENTS::So_hex18PlastType::Instance() or
        actele->ElementType() == DRT::ELEMENTS::So_sh18PlastType::Instance()
       )
    {
      isplasticity_local=true;
      break;
    }

    // Detect EAS --------------------------------------------------------------
    DRT::ELEMENTS::So_base* so_base_ele = dynamic_cast<DRT::ELEMENTS::So_base*>(actele);
    if (so_base_ele!=NULL)
    {
      if (so_base_ele->HaveEAS())
        iseas_local = 1;
    }
    /* additional check for shell8 elements, since these elements are not derived
     * from the So_base class */
    else
    {
      DRT::ELEMENTS::Shell8* shell8_ele = dynamic_cast<DRT::ELEMENTS::Shell8*>(actele);
          if (shell8_ele!=NULL)
            if (shell8_ele->HaveEAS())
              iseas_local = 1;
    }

    // Detect additional pressure dofs -----------------------------------------
    if (actele->ElementType() == DRT::ELEMENTS::So_sh8p8Type::Instance())
      ispressure_local = 1;

    // Detect fbar
    DRT::ELEMENTS::So_hex8fbar* so_hex8fbar_ele =
        dynamic_cast<DRT::ELEMENTS::So_hex8fbar*>(actele);
    if (so_hex8fbar_ele!=NULL)
      isfbar_local = 1;

    // Detect non-additive rotation-vector DOFs --------------------------------
    if (actele->ElementType() == DRT::ELEMENTS::Beam3rType::Instance() or
        actele->ElementType() == DRT::ELEMENTS::Beam3kType::Instance()
       )
    {
      isrotvec_local=true;
      break;
    }
  }

  // plasticity - sum over all processors
  actdis_->Comm().SumAll(&isplasticity_local,&isplasticity_global,1);
  if (isplasticity_global>0)
    eletechs.insert(INPAR::STR::eletech_plasticity);

  // eas - sum over all processors
  actdis_->Comm().SumAll(&iseas_local,&iseas_global,1);
  if (iseas_global>0)
    eletechs.insert(INPAR::STR::eletech_eas);

  // pressure - sum over all processors
  actdis_->Comm().SumAll(&ispressure_local,&ispressure_global,1);
  if (ispressure_global>0)
    eletechs.insert(INPAR::STR::eletech_pressure);

  // fbar - sum over all processors
  actdis_->Comm().SumAll(&isfbar_local,&isfbar_global,1);
  if (isfbar_global>0)
    eletechs.insert(INPAR::STR::eletech_fbar);

  // rotation vector DOFs - sum over all processors
  actdis_->Comm().SumAll(&isrotvec_local,&isrotvec_global,1);
  if (isrotvec_global>0)
    eletechs.insert(INPAR::STR::eletech_rotvec);

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::SetParams(
    Teuchos::ParameterList& ioflags,
    Teuchos::ParameterList& xparams,
    Teuchos::ParameterList& taflags
    )
{
  // get the problem instance and the problem type
  DRT::Problem* problem = DRT::Problem::Instance();
  PROBLEM_TYP probtype = problem->ProblemType();

  // ---------------------------------------------------------------------------
  // show default parameters
  // ---------------------------------------------------------------------------
  if ((actdis_->Comm()).MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(IO::cout, *sdyn_);

  // ---------------------------------------------------------------------------
  // get input parameter lists and copy them,
  // because a few parameters are overwritten
  // ---------------------------------------------------------------------------
  // nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> snox =
      Teuchos::rcp(new Teuchos::ParameterList(problem->StructuralNoxParams()));
  Teuchos::ParameterList& nox = xparams.sublist("NOX");
  nox = *snox;
  // add extra parameters (a kind of work-around)
    xparams.set<FILE*>("err file", problem->ErrorFile()->Handle());
  // Parameter to determine if MLMC is on/off
  Teuchos::RCP<Teuchos::ParameterList> mlmcp
      = Teuchos::rcp(new Teuchos::ParameterList (problem->MultiLevelMonteCarloParams()));
  // Needed for reduced restart output
  xparams.set<int>("REDUCED_OUTPUT",Teuchos::getIntegralValue<int>((*mlmcp),"REDUCED_OUTPUT"));

  sdyn_->set<double>("TIMESTEP", prbdyn_->get<double>("TIMESTEP"));

  // overrule certain parameters
  sdyn_->set<int>("NUMSTEP", prbdyn_->get<int>("NUMSTEP"));
  sdyn_->set<int>("RESTARTEVRY", prbdyn_->get<int>("RESTARTEVRY"));
  sdyn_->set<int>("RESULTSEVRY", prbdyn_->get<int>("RESULTSEVRY"));

  // Check if for chosen Rayleigh damping the regarding parameters are given explicitly in the .dat file
  if (DRT::INPUT::IntegralValue<INPAR::STR::DampKind>(*sdyn_,"DAMPING") == INPAR::STR::damp_rayleigh)
  {
    if (sdyn_->get<double>("K_DAMP") < 0.0)
    {
      dserror("Rayleigh damping parameter K_DAMP not explicitly given.");
    }
    if (sdyn_->get<double>("M_DAMP") < 0.0)
    {
      dserror("Rayleigh damping parameter M_DAMP not explicitly given.");
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  /* Overwrite certain parameters in STRUCTURAL DYNAMIC/TIMADAPTIVITY by those from
   * FSI DYNAMIC/TIMEADAPTIVITY
   *
   * In case, that the structure field is part of an FSI simulation with time step
   * size adaptivity based on structure field error estimation, we have to provide
   * the following algorithmic control parameters:
   *
   * - ADAPTSTEPMAX
   * - STEPSIZEMAX
   * - STEPSIZEMIN
   * - SIZERATIOMAX
   * - SIZERATIOMIN
   * - SIZERATIOSCALE
   *
   * They are specified by the FSI algorithm since they have to be the same for
   * the structure and fluid field. Hence, we overwrite the corresponding
   * parameters in the structural parameter list in order to avoid redundant
   * parameter specification in the input file.
   *
   * Note: This is really ugly, but currently the only way to avoid that the user
   * has to specify these parameters twice in the input file.
   *
   * ToDO: Find something nicer here!
   *
   * \author mayr.mt \date 12/2013
   */
  // ---------------------------------------------------------------------------
  if (probtype == prb_fsi or probtype == prb_fsi_redmodels)
  {
    const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
    const Teuchos::ParameterList& fsiada = fsidyn.sublist("TIMEADAPTIVITY");
    if (DRT::INPUT::IntegralValue<bool>(fsiada,"TIMEADAPTON"))
    {
      // overrule time step size adaptivity control parameters
      if (taflags.get<std::string>("KIND") != "NONE")
      {
        taflags.set<int>("ADAPTSTEPMAX", fsiada.get<int>("ADAPTSTEPMAX"));
        taflags.set<double>("STEPSIZEMAX", fsiada.get<double>("DTMAX"));
        taflags.set<double>("STEPSIZEMIN", fsiada.get<double>("DTMIN"));
        taflags.set<double>("SIZERATIOMAX", fsiada.get<double>("SIZERATIOMAX"));
        taflags.set<double>("SIZERATIOMIN", fsiada.get<double>("SIZERATIOMIN"));
        taflags.set<double>("SIZERATIOSCALE", fsiada.get<double>("SAFETYFACTOR"));

        if (actdis_->Comm().MyPID() == 0)
        {
          IO::cout << "*** Due to FSI time step size adaptivity with structure based error estimation,\n"
                      "algorithmic control parameters in STRUCTURAL DYNAMIC/TIMEADAPTIVITY have been\n"
                      "overwritten by those from FSI DYNAMIC/TIMEADAPTIVITY."
                   << IO::endl << IO::endl;
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::SetStructureWrapper(
    const Teuchos::ParameterList& ioflags,
    const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams,
    const Teuchos::ParameterList& taflags,
    Teuchos::RCP<STR::TIMINT::Base> ti_strategy)
{
  // create a adaptive wrapper
  CreateAdaptiveWrapper(ioflags,sdyn,xparams,taflags,ti_strategy);

  // if no adaptive wrapper was found, we try to create a standard one
  if (str_wrapper_.is_null())
    CreateWrapper(ti_strategy);

  if (str_wrapper_.is_null())
    dserror("No proper time integration found!");

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::CreateAdaptiveWrapper(
    const Teuchos::ParameterList& ioflags,
    const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams,
    const Teuchos::ParameterList& taflags,
    Teuchos::RCP<STR::TIMINT::Base> ti_strategy)
{
  // get the problem instance and the problem type
  DRT::Problem* problem = DRT::Problem::Instance();
  PROBLEM_TYP probtype = problem->ProblemType();

  // create auxiliary time integrator, can be seen as a wrapper for ti_strategy
  Teuchos::RCP<STR::TimAda> wrapper_adaptive =
      STR::TIMINT::BuildAdaptiveWrapper(ioflags, sdyn, xparams, taflags, ti_strategy);

  if (wrapper_adaptive.is_null())
    return;

  switch (probtype)
  {
    case prb_structure: // pure structural time adaptivity
    case prb_crack:
    {
      str_wrapper_ = Teuchos::rcp(new StructureTimIntAda(wrapper_adaptive, ti_strategy));
      break;
    }
    case prb_fsi: // structure based time adaptivity within an FSI simulation
    case prb_fsi_redmodels:
    {
      if ((actdis_->Comm()).MyPID()==0)
        IO::cout << "Using StructureNOXCorrectionWrapper()..." << IO::endl;

      Teuchos::RCP<FSIStructureWrapper> fsiwrapperwithadaptivity =
          Teuchos::rcp(new StructureFSITimIntAda(wrapper_adaptive, Teuchos::rcp(new StructureNOXCorrectionWrapper(ti_strategy))));
      str_wrapper_ = fsiwrapperwithadaptivity;
      break;
    }
    default:
    {
      dserror("Adaptive time integration for the structure not implemented for desired problem type.");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithmNew::CreateWrapper(
    Teuchos::RCP<STR::TIMINT::Base> ti_strategy)
{
  // get the problem instance and the problem type
  DRT::Problem* problem = DRT::Problem::Instance();
  PROBLEM_TYP probtype = problem->ProblemType();

  switch(probtype)
  {
    case prb_fsi:
    case prb_fsi_redmodels:
    case prb_fsi_lung:
    case prb_gas_fsi:
    case prb_ac_fsi:
    case prb_biofilm_fsi:
    case prb_thermo_fsi:
    case prb_fsi_xfem:
    {
      const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
      const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");

      // Are there any constraint conditions active?
      const std::set<INPAR::STR::ModelType>& modeltypes = ti_strategy->GetDataSDyn().GetModelTypes();
      if (modeltypes.find(INPAR::STR::model_lag_pen_constraint)!=modeltypes.end())
      {
        if ((actdis_->Comm()).MyPID()==0)
          IO::cout << "Using StructureNOXCorrectionWrapper()..." << IO::endl;

        if (coupling == fsi_iter_constr_monolithicstructuresplit or
            coupling == fsi_iter_constr_monolithicfluidsplit)
          str_wrapper_ = Teuchos::rcp(new FSIStructureWrapper(Teuchos::rcp(new StructureNOXCorrectionWrapper(ti_strategy))));
        else
          str_wrapper_ = Teuchos::rcp(new StructureConstrMerged(Teuchos::rcp(new StructureNOXCorrectionWrapper(ti_strategy))));
      }
      else
      {
        if (coupling == fsi_iter_lung_monolithicstructuresplit or
            coupling == fsi_iter_lung_monolithicfluidsplit)
        {
          if ((actdis_->Comm()).MyPID()==0)
            IO::cout << "Using StructureNOXCorrectionWrapper()..." << IO::endl;
          str_wrapper_ = Teuchos::rcp(new StructureLung(Teuchos::rcp(new StructureNOXCorrectionWrapper(ti_strategy))));
        }
        else
        {
          if(problem->Restart())
            ti_strategy->Setup();
          str_wrapper_ = Teuchos::rcp(new FSIStructureWrapper(ti_strategy)); // case of partitioned fsi
        }
      }
      break;
    }
    case prb_immersed_fsi:
    case prb_immersed_ale_fsi:
    {
      if(problem->Restart())
        ti_strategy->Setup();
      str_wrapper_ = Teuchos::rcp(new FSIStructureWrapperImmersed(ti_strategy));
      break;
    }
    case prb_fsi_crack:
      str_wrapper_ = Teuchos::rcp(new FSICrackingStructure(Teuchos::rcp(new FSIStructureWrapper(Teuchos::rcp(new StructureNOXCorrectionWrapper(ti_strategy))))));
      break;
    case prb_redairways_tissue:
      str_wrapper_ = Teuchos::rcp(new StructureRedAirway(ti_strategy));
      break;
    case prb_poroelast:
    case prb_poroscatra:
    case prb_fpsi:
    case prb_fps3i:
    case prb_fpsi_xfem:
    case prb_immersed_cell:
    {
      const Teuchos::ParameterList& porodyn = problem->PoroelastDynamicParams();
      const INPAR::POROELAST::SolutionSchemeOverFields coupling =
            DRT::INPUT::IntegralValue<INPAR::POROELAST::SolutionSchemeOverFields>(porodyn, "COUPALGO");
      // Are there any constraint conditions active?
      const std::set<INPAR::STR::ModelType>& modeltypes = ti_strategy->GetDataSDyn().GetModelTypes();
      if (modeltypes.find(INPAR::STR::model_lag_pen_constraint)!=modeltypes.end())
      {
        if (   coupling == INPAR::POROELAST::Monolithic_structuresplit
            or coupling == INPAR::POROELAST::Monolithic_fluidsplit
            or coupling == INPAR::POROELAST::Monolithic_nopenetrationsplit
            )
          str_wrapper_ = Teuchos::rcp(new FPSIStructureWrapper(ti_strategy));
        else
          str_wrapper_ = Teuchos::rcp(new StructureConstrMerged(ti_strategy));
      }
      else
      {
          str_wrapper_ = Teuchos::rcp(new FPSIStructureWrapper(ti_strategy));
      }
      break;
    }
    case prb_struct_ale:
      str_wrapper_ = Teuchos::rcp(new FSIStructureWrapper(ti_strategy));
      break;
    case prb_invana:
      str_wrapper_ = (Teuchos::rcp(new StructureInvana(ti_strategy)));
      break;
    default:
      /// wrap time loop for pure structure problems
      str_wrapper_ = (Teuchos::rcp(new StructureTimeLoop(ti_strategy)));
      break;
  }

  return;
}


