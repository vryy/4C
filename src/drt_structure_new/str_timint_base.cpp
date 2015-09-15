/*
 * str_timint_base.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */


#include "str_timint_base.H"
#include "str_model_evaluator_factory.H"
#include "str_model_evaluator_generic.H"

#include "../drt_lib/drt_utils_timintmstep.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Epetra_Map.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataIO::BaseDataIO()
    : isinit_(false),
      issetup_(false),
      output_(Teuchos::null),
      energyfile_(Teuchos::null),
      errfile_(NULL),
      gmsh_out_(false),
      printlogo_(false),
      printerrfile_(false),
      outputeveryiter_(false),
      writesurfactant_(false),
      writestate_(false),
      writevelacc_(false),
      printscreen_(-1),
      oei_filecounter_(-1),
      outputcounter_(-1),
      writerestartevery_(-1),
      writereducedrestart_(-1),
      writeresultsevery_(-1),
      writeenergyevery_(-1),
      kinergy_(-1.0),
      intergy_(-1.0),
      extergy_(-1.0),
      writestress_(INPAR::STR::stress_none),
      writecouplstress_(INPAR::STR::stress_none),
      writestrain_(INPAR::STR::strain_none),
      writeplstrain_(INPAR::STR::strain_none)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::Init(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<IO::DiscretizationWriter> output
    )
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // ----------------------------------------------------------
  // initialize the printing and output parameters
  // ----------------------------------------------------------
  {
    output_ = output;
    printscreen_ = ioparams.get<int>("STDOUTEVRY");
    printlogo_ = (printscreen_>0 ? true : false);
    errfile_ = xparams.get<FILE*>("err file");
    printerrfile_ = (true and errfile_);
    printiter_ = true;
    outputeveryiter_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"OUTPUT_EVERY_ITER");
    oei_filecounter_ = ioparams.get<int>("OEI_FILE_COUNTER");
    writerestartevery_ = sdynparams.get<int>("RESTARTEVRY");
    writereducedrestart_ = xparams.get<int>("REDUCED_OUTPUT");
    writestate_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_DISP");
    writevelacc_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_VEL_ACC");
    writeresultsevery_ = sdynparams.get<int>("RESULTSEVRY");
    writestress_ = DRT::INPUT::IntegralValue<INPAR::STR::StressType>(ioparams,"STRUCT_STRESS");
    writecouplstress_ = DRT::INPUT::IntegralValue<INPAR::STR::StressType>(ioparams,"STRUCT_COUPLING_STRESS");
    writestrain_ = DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(ioparams,"STRUCT_STRAIN");
    writeplstrain_ = DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(ioparams,"STRUCT_PLASTIC_STRAIN");
    writeenergyevery_ = sdynparams.get<int>("RESEVRYERGY");
    writesurfactant_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_SURFACTANT");
  }

  isinit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  issetup_ = true;

  // Good bye
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataSDyn::BaseDataSDyn()
    : isinit_(false),
      issetup_(false),
      timemax_(-1.0),
      stepmax_(-1),
      timer_(Teuchos::null),
      damptype_(INPAR::STR::damp_none),
      dampk_(-1.0),
      dampm_(-1.0),
      masslintype_(INPAR::STR::ml_none),
      modeltypes_(Teuchos::null),
      dyntype_(INPAR::STR::dyna_statics),
      prestresstype_(INPAR::STR::prestress_none),
      nlnsolvertype_(INPAR::STR::soltech_vague),
      noxparams_(Teuchos::null),
      linsolvers_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataSDyn::Init(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::ParameterList& sDynParams,
    const Teuchos::ParameterList& xparams,
    const Teuchos::RCP<std::vector<const enum INPAR::STR::ModelType> > modeltypes,
    const Teuchos::RCP<std::map<const enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > > linsolvers
    )
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // ----------------------------------------------------------
  // initialize general variables
  // ----------------------------------------------------------
  {
    timemax_ = sDynParams.get<double>("MAXTIME");
    stepmax_ = sDynParams.get<int>("NUMSTEP");

    timer_ = Teuchos::rcp(new Epetra_Time(discret->Comm()));

    dyntype_ =
        DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sDynParams, "DYNAMICTYP");
  }
  // ----------------------------------------------------------
  // initialize the damping control parameters
  // ----------------------------------------------------------
  {
    damptype_ = DRT::INPUT::IntegralValue<INPAR::STR::DampKind>(sDynParams,"DAMPING");
    dampk_ = sDynParams.get<double>("K_DAMP");
    dampm_ = sDynParams.get<double>("M_DAMP");
  }
  // ----------------------------------------------------------
  // initialize the mass and inertia control parameters
  // ----------------------------------------------------------
  {
    masslintype_ = DRT::INPUT::IntegralValue<INPAR::STR::MassLin>(sDynParams,"MASSLIN");
  }
  // ----------------------------------------------------------
  // initialize model evaluator control parameters
  // ----------------------------------------------------------
  {
    modeltypes_ = modeltypes;
  }
  // ----------------------------------------------------------
  // initialize implicit variables
  // ----------------------------------------------------------
  {
    prestresstype_ =
          DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sDynParams,"PRESTRESS");
    nlnsolvertype_ =
        DRT::INPUT::IntegralValue<INPAR::STR::NonlinSolTech>(sDynParams,"NLNSOL");
    noxparams_ = Teuchos::rcp(new Teuchos::ParameterList(xparams.sublist("NOX")));
  }
  // ----------------------------------------------------------
  // initialize linear solver variables
  // ----------------------------------------------------------
  {
    linsolvers_ = linsolvers;
  }

  isinit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataSDyn::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  issetup_ = true;

  // Good bye
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState::BaseDataGlobalState()
    : isinit_(false),
      issetup_(false),
      datasdyn_(Teuchos::null),
      discret_(Teuchos::null),
      comm_(Teuchos::null),
      myRank_(-1),
      timenp_(0.0),
      timen_(Teuchos::null),
      dt_(Teuchos::null),
      stepn_(0),
      stepnp_(0),
      dis_(Teuchos::null),
      vel_(Teuchos::null),
      acc_(Teuchos::null),
      disnp_(Teuchos::null),
      velnp_(Teuchos::null),
      accnp_(Teuchos::null),
      stiff_(Teuchos::null),
      mass_(Teuchos::null),
      damp_(Teuchos::null),
      timer_(Teuchos::null),
      dtsolve_(0.0),
      dtele_(0.0)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataGlobalState::Init(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<const BaseDataSDyn> datasdyn
    )
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // ----------------------------------------------------------
  // initialize a const pointer to the sDynData container
  // ----------------------------------------------------------
  datasdyn_ = datasdyn;

  // ----------------------------------------------------------
  // initialize general purpose algorithm members
  // ----------------------------------------------------------
  {
    discret_ = discret;
    comm_ = Teuchos::rcpFromRef(discret_->Comm());
    myRank_  = comm_->MyPID();
  }

  // end of initialization
  isinit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataGlobalState::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  // --------------------------------------
  // setup state vectors
  // --------------------------------------
  // displacements D_{n}
  dis_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  // velocities V_{n}
  vel_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  // accelerations A_{n}
  acc_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));

  // displacements D_{n+1} at t_{n+1}
  disnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  // velocities V_{n+1} at t_{n+1}
  velnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  // accelerations A_{n+1} at t_{n+1}
  accnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  // --------------------------------------
  // setup sparse operators
  // --------------------------------------
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, true, true));
  mass_  = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, true, true));
  if (datasdyn_->GetDampingType() != INPAR::STR::damp_none)
  {
    if (datasdyn_->GetMassLinType() == INPAR::STR::ml_none)
    {
      damp_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, true, true));
    }
    else
    {
      //Since our element evaluate routine is only designed for two input matrices
      //(stiffness and damping or stiffness and mass) its not possible, to have nonlinear
      //inertia forces AND material damping.
      dserror("So far it is not possible to model nonlinear inertia forces and damping!");
    }
  }

  issetup_ = true;

  // Good bye
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::TIMINT::BaseDataGlobalState::DofRowMap() const
{
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::TIMINT::BaseDataGlobalState::DofRowMap(unsigned nds) const
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap(nds);
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map* STR::TIMINT::BaseDataGlobalState::DofRowMapView() const
{
  return discret_->DofRowMap();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Base::Base()
    : isinit_(false),
      issetup_(false),
      dataio_(Teuchos::null),
      datasdyn_(Teuchos::null),
      dataglobalstate_(Teuchos::null),
      modelevaluators_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate
    )
{
  // ---------------------------------------------
  // We need to call Setup() after Init()
  // ---------------------------------------------
  issetup_ = false;

  // ---------------------------------------------
  // initilize the data container ptrs
  // ---------------------------------------------
  dataio_ = dataio;
  datasdyn_ = datasdyn;
  dataglobalstate_ = dataglobalstate;

  // ---------------------------------------------
  // build model evaluator
  // ---------------------------------------------
  modelevaluators_ = STR::MODELEVALUATOR::BuildModelEvaluators(DataSDyn().GetModelType());

  std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >::iterator me_iter;
  for (me_iter=ModelEvaluatorsPtr()->begin();me_iter!=ModelEvaluatorsPtr()->end();++me_iter)
  {
    me_iter->second->Init();
    me_iter->second>Setup();
  }

  // ---------------------------------------------
  // set isInit flag
  // ---------------------------------------------
  isinit_ = true;

  // good bye
  return;
}


/*----------------------------------------------------------------------*/
/* Read and set restart values */
void STR::TIMINT::Base::ReadRestart
(
  const int step
)
{
  IO::DiscretizationReader reader(DataGlobalStatePtr()->GetDiscret(), step);
  if (step != reader.ReadInt("step"))
    dserror("Time step on file not equal to given step");

//  step_ = step;
//  stepn_ = step_ + 1;
//  time_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
//  timen_ = (*time_)[0] + (*dt_)[0];


  // TODO: restart for model--evaluators...

//  ReadRestartState();
//
//  ReadRestartConstraint();
//  ReadRestartWindkessel();
//  ReadRestartContactMeshtying();
//  ReadRestartBeamContact();
//  ReadRestartStatMech();
//  ReadRestartSurfstress();
//  ReadRestartMultiScale();
//  ReadRestartCrack();
//
//  ReadRestartForce();

}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::TIMINT::Base::DofRowMap()
{
  CheckInitSetup();
  const Epetra_Map* dofrowmap = DataGlobalStatePtr()->GetDiscret()->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::Base::SystemMatrix()
{
  CheckInitSetup();
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(DataGlobalState().GetMutableStiffMatrix());
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> STR::TIMINT::Base::BlockSystemMatrix()
{
  CheckInitSetup();
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(DataGlobalState().GetMutableStiffMatrix());
}
