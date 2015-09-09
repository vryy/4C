/*
 * str_timint_base.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */


#include "str_timint_base.H"

#include "../drt_lib/drt_utils_timintmstep.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Time.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState::BaseDataGlobalState()
    : isInit_(false),
      isSetup_(false),
      dataSDyn_(Teuchos::null),
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
    const Teuchos::RCP<const BaseDataSDyn> dataSDyn
    )
{
  // We have to call Setup() after Init()
  isSetup_ = false;

  // ----------------------------------------------------------
  // initialize a const pointer to the sDynData container
  // ----------------------------------------------------------
  dataSDyn_ = dataSDyn;

  // ----------------------------------------------------------
  // initialize general purpose algorithm members
  // ----------------------------------------------------------
  {
    discret_ = discret;
    comm_ = Teuchos::rcpFromRef(discret_->Comm());
    myRank_  = comm_->MyPID();
  }

  // end of initialization
  isInit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataGlobalState::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  // displacements D_{n}
  dis_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, discret_->DofRowMap(), true));
  // velocities V_{n}
  vel_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, discret_->DofRowMap(), true));
  // accelerations A_{n}
  acc_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, discret_->DofRowMap(), true));

  // displacements D_{n+1} at t_{n+1}
  disnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  // velocities V_{n+1} at t_{n+1}
  velnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  // accelerations A_{n+1} at t_{n+1}
  accnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  // create empty matrices
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*discret_->DofRowMap(), 81, true, true));
  mass_  = Teuchos::rcp(new LINALG::SparseMatrix(*discret_->DofRowMap(), 81, true, true));
  if (dataSDyn_->GetDampingType() != INPAR::STR::damp_none)
  {
    if (dataSDyn_->GetMassLinType() == INPAR::STR::ml_none)
    {
      damp_ = Teuchos::rcp(new LINALG::SparseMatrix(*discret_->DofRowMap(), 81, true, true));
    }
    else
    {
      //Since our element evaluate routine is only designed for two input matrices
      //(stiffness and damping or stiffness and mass) its not possible, to have nonlinear
      //inertia forces AND material damping.
      dserror("So far it is not possible to model nonlinear inertia forces and damping!");
    }
  }

  isSetup_ = true;

  // Good bye
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataIO::BaseDataIO()
    : isInit_(false),
      isSetup_(false),
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
void STR::TIMINT::BaseDataIO::Init(const Teuchos::ParameterList& IOParams,
    const Teuchos::ParameterList& sDynParams,
    const Teuchos::ParameterList& xParams,
    Teuchos::RCP<IO::DiscretizationWriter> output
    )
{
  // We have to call Setup() after Init()
  isSetup_ = false;

  // ----------------------------------------------------------
  // initialize the printing and output parameters
  // ----------------------------------------------------------
  {
    output_ = output;
    printscreen_ = IOParams.get<int>("STDOUTEVRY");
    printlogo_ = (printscreen_>0 ? true : false);
    errfile_ = xParams.get<FILE*>("err file");
    printerrfile_ = (true and errfile_);
    printiter_ = true;
    outputeveryiter_ = (bool) DRT::INPUT::IntegralValue<int>(IOParams,"OUTPUT_EVERY_ITER");
    oei_filecounter_ = IOParams.get<int>("OEI_FILE_COUNTER");
    writerestartevery_ = sDynParams.get<int>("RESTARTEVRY");
    writereducedrestart_ = xParams.get<int>("REDUCED_OUTPUT");
    writestate_ = (bool) DRT::INPUT::IntegralValue<int>(IOParams,"STRUCT_DISP");
    writevelacc_ = (bool) DRT::INPUT::IntegralValue<int>(IOParams,"STRUCT_VEL_ACC");
    writeresultsevery_ = sDynParams.get<int>("RESULTSEVRY");
    writestress_ = DRT::INPUT::IntegralValue<INPAR::STR::StressType>(IOParams,"STRUCT_STRESS");
    writecouplstress_ = DRT::INPUT::IntegralValue<INPAR::STR::StressType>(IOParams,"STRUCT_COUPLING_STRESS");
    writestrain_ = DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(IOParams,"STRUCT_STRAIN");
    writeplstrain_ = DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(IOParams,"STRUCT_PLASTIC_STRAIN");
    writeenergyevery_ = sDynParams.get<int>("RESEVRYERGY");
    writesurfactant_ = (bool) DRT::INPUT::IntegralValue<int>(IOParams,"STRUCT_SURFACTANT");
  }

  isInit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  isSetup_ = true;

  // Good bye
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataSDyn::BaseDataSDyn()
    : isInit_(false),
      isSetup_(false),
      timemax_(-1.0),
      stepmax_(-1),
      timer_(Teuchos::null),
      damptype_(INPAR::STR::damp_none),
      dampk_(-1.0),
      dampm_(-1.0),
      masslintype_(INPAR::STR::ml_none),
      modeltypes_(Teuchos::null),
      dynType_(INPAR::STR::dyna_statics),
      preStressType_(INPAR::STR::prestress_none),
      nlnSolverType_(INPAR::STR::soltech_vague)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataSDyn::Init(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::ParameterList& sDynParams,
    const Teuchos::RCP<std::vector<const enum INPAR::STR::ModelType> > modeltypes
    )
{
  // We have to call Setup() after Init()
  isSetup_ = false;

  // ----------------------------------------------------------
  // initialize general variables
  // ----------------------------------------------------------
  {
    timemax_ = sDynParams.get<double>("MAXTIME");
    stepmax_ = sDynParams.get<int>("NUMSTEP");

    timer_ = Teuchos::rcp(new Epetra_Time(discret->Comm()));

    dynType_ =
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
    masslintype_ = DRT::INPUT::IntegralValue<INPAR::STR::DampKind>(sDynParams,"MASSLIN");
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
    preStressType_ =
          DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sDynParams,"PRESTRESS");
    nlnSolverType_ =
        DRT::INPUT::IntegralValue<INPAR::STR::NonlinSolTech>(sDynParams,"NLNSOL");
  }

  isInit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataSDyn::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  isSetup_ = true;

  // Good bye
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Base::Base() :
  isInit_(false),
  isSetup_(false),
  dataIO_(Teuchos::null),
  dataSDyn_(Teuchos::null),
  dataGlobalState_(Teuchos::null),
  linsolvers_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate,
    const Teuchos::RCP<std::map<const enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > > linsolvers
    )
{
  // We need to call Setup() after Init()
  isSetup_ = false;

  // initilize the data container ptrs
  dataIO_ = dataio;
  dataSDyn_ = datasdyn;
  dataGlobalState_ = dataglobalstate;

  // initialize the linear solvers ptr
  linsolvers_ = linsolvers;

  // set isInit flag
  isInit_ = true;

  // good bye
  return;
}
