/*
 * str_timint_base.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */


#include "str_timint_base.H"
#include "str_model_evaluator_factory.H"
#include "str_model_evaluator_generic.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_timestepping/timintmstep.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"

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
      printiter_(false),
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

  // ---------------------------------------------------------------------------
  // initialize the printing and output parameters
  // ---------------------------------------------------------------------------
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
      eletechs_(Teuchos::null),
      dyntype_(INPAR::STR::dyna_statics),
      itermin_(-1),
      itermax_(-1),
      prestresstype_(INPAR::STR::prestress_none),
      predtype_(INPAR::STR::pred_vague),
      nlnsolvertype_(INPAR::STR::soltech_vague),
      divergenceaction_(INPAR::STR::divcont_stop),
      noxparams_(Teuchos::null),
      linsolvers_(Teuchos::null),
      normtype_(INPAR::STR::norm_vague),
      tol_disp_incr_(-1.0),
      toltype_disp_incr_(INPAR::STR::convnorm_abs),
      tol_fres_(-1.0),
      toltype_fres_(INPAR::STR::convnorm_abs),
      tol_pres_(-1.0),
      toltype_pres_(INPAR::STR::convnorm_abs),
      tol_inco_(-1.0),
      toltype_inco_(INPAR::STR::convnorm_abs),
      tol_plast_res_(-1.0),
      toltype_plast_res_(INPAR::STR::convnorm_abs),
      tol_plast_incr_(-1.0),
      toltype_plast_incr_(INPAR::STR::convnorm_abs),
      tol_eas_res_(-1.0),
      toltype_eas_res_(INPAR::STR::convnorm_abs),
      tol_eas_incr_(-1.0),
      toltype_eas_incr_(INPAR::STR::convnorm_abs),
      normcombo_disp_pres_(INPAR::STR::bop_and),
      normcombo_fres_inco_(INPAR::STR::bop_and),
      normcombo_fres_eas_res_(INPAR::STR::bop_and),
      normcombo_disp_eas_incr_(INPAR::STR::bop_and),
      normcombo_fres_plast_res_(INPAR::STR::bop_and),
      normcombo_disp_plast_incr_(INPAR::STR::bop_and),
      normcombo_fres_disp_(INPAR::STR::bop_and),
      toltype_constr_res_(INPAR::STR::convnorm_abs),
      tol_constr_res_(-1.0),
      toltype_windk_res_(INPAR::STR::convnorm_abs),
      tol_windk_res_(-1.0),
      toltype_windk_incr_(INPAR::STR::convnorm_abs),
      tol_windk_incr_(-1.0),
      toltype_contact_res_(INPAR::STR::convnorm_abs),
      tol_contact_res_(-1.0),
      toltype_contact_lm_incr_(INPAR::STR::convnorm_abs),
      tol_contact_lm_incr_(-1.0),
      normcombo_fres_contact_res_(INPAR::STR::bop_and),
      normcombo_disp_contact_lm_incr_(INPAR::STR::bop_and)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataSDyn::Init(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams,
    const Teuchos::RCP<std::set<enum INPAR::STR::ModelType> > modeltypes,
    const Teuchos::RCP<std::set<enum INPAR::STR::EleTech> > eletechs,
    const Teuchos::RCP<std::map<enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > > linsolvers
    )
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // ---------------------------------------------------------------------------
  // initialize general variables
  // ---------------------------------------------------------------------------
  {
    timemax_ = sdynparams.get<double>("MAXTIME");
    stepmax_ = sdynparams.get<int>("NUMSTEP");

    timer_ = Teuchos::rcp(new Epetra_Time(discret->Comm()));

    dyntype_ =
        DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP");
  }
  // ---------------------------------------------------------------------------
  // initialize the damping control parameters
  // ---------------------------------------------------------------------------
  {
    damptype_ = DRT::INPUT::IntegralValue<INPAR::STR::DampKind>(sdynparams,"DAMPING");
    dampk_ = sdynparams.get<double>("K_DAMP");
    dampm_ = sdynparams.get<double>("M_DAMP");
  }
  // ---------------------------------------------------------------------------
  // initialize the mass and inertia control parameters
  // ---------------------------------------------------------------------------
  {
    masslintype_ = DRT::INPUT::IntegralValue<INPAR::STR::MassLin>(sdynparams,"MASSLIN");
  }
  // ---------------------------------------------------------------------------
  // initialize model evaluator control parameters
  // ---------------------------------------------------------------------------
  {
    modeltypes_ = modeltypes;
    eletechs_   = eletechs;
  }
  // ---------------------------------------------------------------------------
  // initialize implicit variables
  // ---------------------------------------------------------------------------
  {
    itermin_ = sdynparams.get<int>("MAXITER");
    itermax_ = sdynparams.get<int>("MINITER");
    prestresstype_ =
          DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdynparams,"PRESTRESS");
    nlnsolvertype_ =
        DRT::INPUT::IntegralValue<INPAR::STR::NonlinSolTech>(sdynparams,"NLNSOL");
    noxparams_ = Teuchos::rcp(new Teuchos::ParameterList(xparams.sublist("NOX")));
  }
  // ---------------------------------------------------------------------------
  // initialize linear solver variables
  // ---------------------------------------------------------------------------
  {
    linsolvers_ = linsolvers;
  }
  // ---------------------------------------------------------------------------
  // initialize the status test control parameters
  // ---------------------------------------------------------------------------
  {
    normtype_ = DRT::INPUT::IntegralValue<INPAR::STR::VectorNorm>(sdynparams,"ITERNORM");

    // -------------------------------------------------------------------------
    // primary variables
    // -------------------------------------------------------------------------
    tol_disp_incr_ = sdynparams.get<double>("TOLDISP");
    toltype_disp_incr_ = DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(
        sdynparams,"NORM_DISP");

    tol_fres_ = sdynparams.get<double>("TOLRES");
    toltype_fres_ = DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(
        sdynparams,"NORM_RESF");

    tol_pres_ = sdynparams.get<double>("TOLPRE");
    toltype_pres_ = INPAR::STR::convnorm_abs;

    tol_inco_ = sdynparams.get<double>("TOLINCO");
    toltype_inco_ = INPAR::STR::convnorm_abs;

    tol_plast_res_ = DRT::Problem::Instance()->
        SemiSmoothPlastParams().get<double>("TOLPLASTCONSTR");
    toltype_plast_res_ = INPAR::STR::convnorm_abs;

    tol_plast_incr_ = DRT::Problem::Instance()->
        SemiSmoothPlastParams().get<double>("TOLDELTALP");
    toltype_plast_incr_ = INPAR::STR::convnorm_abs;

    tol_eas_res_ = DRT::Problem::Instance()->
        SemiSmoothPlastParams().get<double>("TOLEASRES");
    toltype_eas_res_ = INPAR::STR::convnorm_abs;

    tol_eas_incr_ = DRT::Problem::Instance()->
        SemiSmoothPlastParams().get<double>("TOLEASINCR");
    toltype_eas_incr_ = INPAR::STR::convnorm_abs;

    normcombo_disp_pres_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        sdynparams,"NORMCOMBI_DISPPRES");
    normcombo_fres_inco_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        sdynparams,"NORMCOMBI_RESFINCO");
    normcombo_fres_plast_res_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        DRT::Problem::Instance()->SemiSmoothPlastParams(),"NORMCOMBI_RESFPLASTCONSTR");
    normcombo_disp_plast_incr_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        DRT::Problem::Instance()->SemiSmoothPlastParams(),"NORMCOMBI_DISPPLASTINCR");
    normcombo_fres_eas_res_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        DRT::Problem::Instance()->SemiSmoothPlastParams(),"NORMCOMBI_EASRES");
    normcombo_disp_eas_incr_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        DRT::Problem::Instance()->SemiSmoothPlastParams(),"NORMCOMBI_EASINCR");
    normcombo_fres_disp_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        sdynparams,"NORMCOMBI_RESFDISP");

    // -------------------------------------------------------------------------
    // constraint variables
    // -------------------------------------------------------------------------
    tol_constr_res_ = sdynparams.get<double>("TOLCONSTR");
    toltype_constr_res_ = INPAR::STR::convnorm_abs;

    tol_windk_res_ = DRT::Problem::Instance()->
        WindkesselStructuralParams().get<double>("TOLWINDKESSEL");
    toltype_windk_res_ = INPAR::STR::convnorm_abs;

    tol_windk_incr_ = DRT::Problem::Instance()->
        WindkesselStructuralParams().get<double>("TOLWINDKESSELDOFINCR");
    toltype_windk_incr_ = INPAR::STR::convnorm_abs;

    tol_contact_res_ = DRT::Problem::Instance()->
        ContactDynamicParams().get<double>("TOLCONTCONSTR");
    toltype_contact_res_ = INPAR::STR::convnorm_abs;

    tol_contact_lm_incr_ = DRT::Problem::Instance()->
        ContactDynamicParams().get<double>("TOLLAGR");
    toltype_contact_lm_incr_ = INPAR::STR::convnorm_abs;

    normcombo_fres_contact_res_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        DRT::Problem::Instance()->ContactDynamicParams(),"NORMCOMBI_RESFCONTCONSTR");
    normcombo_disp_contact_lm_incr_ = DRT::INPUT::IntegralValue<INPAR::STR::BinaryOp>(
        DRT::Problem::Instance()->ContactDynamicParams(),"NORMCOMBI_DISPLAGR");
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
double STR::TIMINT::BaseDataSDyn::GetResTolerance(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  switch (qtype)
  {
    case NOX::NLN::StatusTest::quantity_structure:
      return tol_fres_;
      break;
    case NOX::NLN::StatusTest::quantity_contact:
    case NOX::NLN::StatusTest::quantity_meshtying:
      return tol_contact_res_;
      break;
    case NOX::NLN::StatusTest::quantity_windkessel:
      return tol_windk_res_;
      break;
    case NOX::NLN::StatusTest::quantity_lag_pen_constraint:
      return tol_constr_res_;
      break;
    case NOX::NLN::StatusTest::quantity_plasticity:
      return tol_plast_res_;
      break;
    case NOX::NLN::StatusTest::quantity_pressure:
      return tol_inco_;
      break;
    case NOX::NLN::StatusTest::quantity_eas:
      return tol_eas_res_;
      break;
    default:
      dserror("There is no residual tolerance for the given quantity type! "
          "(quantity: %s)",NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      break;
  }

  return -1.0;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::BaseDataSDyn::GetIncrTolerance(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  switch (qtype)
  {
    case NOX::NLN::StatusTest::quantity_structure:
      return tol_disp_incr_;
      break;
    case NOX::NLN::StatusTest::quantity_contact:
    case NOX::NLN::StatusTest::quantity_meshtying:
      return tol_contact_lm_incr_;
      break;
    case NOX::NLN::StatusTest::quantity_windkessel:
      return tol_windk_incr_;
      break;
    case NOX::NLN::StatusTest::quantity_plasticity:
      return tol_plast_incr_;
      break;
    case NOX::NLN::StatusTest::quantity_pressure:
      return tol_pres_;
      break;
    case NOX::NLN::StatusTest::quantity_eas:
      return tol_eas_incr_;
      break;
    default:
      dserror("There is no increment tolerance for the given quantity type! "
          "(quantity: %s)",NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      break;
  }

  return -1.0;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::ConvNorm STR::TIMINT::BaseDataSDyn::GetResToleranceType(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  switch (qtype)
  {
    case NOX::NLN::StatusTest::quantity_structure:
      return toltype_fres_;
      break;
    case NOX::NLN::StatusTest::quantity_contact:
    case NOX::NLN::StatusTest::quantity_meshtying:
      return toltype_contact_res_;
      break;
    case NOX::NLN::StatusTest::quantity_windkessel:
      return toltype_windk_res_;
      break;
    case NOX::NLN::StatusTest::quantity_lag_pen_constraint:
      return toltype_constr_res_;
      break;
    case NOX::NLN::StatusTest::quantity_plasticity:
      return toltype_plast_res_;
      break;
    case NOX::NLN::StatusTest::quantity_pressure:
      return toltype_inco_;
      break;
    case NOX::NLN::StatusTest::quantity_eas:
      return toltype_eas_res_;
      break;
    default:
      dserror("There is no residual tolerance type for the given quantity type! "
          "(quantity: %s)",NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      break;
  }

  return INPAR::STR::convnorm_abs;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::ConvNorm STR::TIMINT::BaseDataSDyn::GetIncrToleranceType(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  switch (qtype)
  {
    case NOX::NLN::StatusTest::quantity_structure:
      return toltype_disp_incr_;
      break;
    case NOX::NLN::StatusTest::quantity_contact:
    case NOX::NLN::StatusTest::quantity_meshtying:
      return toltype_contact_lm_incr_;
      break;
    case NOX::NLN::StatusTest::quantity_windkessel:
      return toltype_windk_incr_;
      break;
    case NOX::NLN::StatusTest::quantity_plasticity:
      return toltype_plast_incr_;
      break;
    case NOX::NLN::StatusTest::quantity_pressure:
      return toltype_pres_;
      break;
    case NOX::NLN::StatusTest::quantity_eas:
      return toltype_eas_incr_;
      break;
    default:
      dserror("There is no increment tolerance type for the given quantity type! "
          "(quantity: %s)",NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      break;
  }

  return INPAR::STR::convnorm_abs;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::BinaryOp STR::TIMINT::BaseDataSDyn::GetResComboType(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  return GetResComboType(NOX::NLN::StatusTest::quantity_structure,qtype);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::BinaryOp STR::TIMINT::BaseDataSDyn::GetResComboType(
    const enum NOX::NLN::StatusTest::QuantityType& qtype_1,
    const enum NOX::NLN::StatusTest::QuantityType& qtype_2) const
{
  // combination: STRUCTURE <--> PRESSURE
  if ((qtype_1==NOX::NLN::StatusTest::quantity_structure and
       qtype_2==NOX::NLN::StatusTest::quantity_pressure) or
      (qtype_1==NOX::NLN::StatusTest::quantity_pressure and
       qtype_2==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_fres_inco_;
  // combination: STRUCTURE <--> EAS
  else if ((qtype_1==NOX::NLN::StatusTest::quantity_structure and
        qtype_2==NOX::NLN::StatusTest::quantity_eas) or
       (qtype_1==NOX::NLN::StatusTest::quantity_eas and
        qtype_2==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_fres_eas_res_;
  // combination: STRUCTURE <--> PLASTICITY
  else if ((qtype_1==NOX::NLN::StatusTest::quantity_structure and
        qtype_2==NOX::NLN::StatusTest::quantity_plasticity) or
       (qtype_1==NOX::NLN::StatusTest::quantity_plasticity and
        qtype_2==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_fres_plast_res_;
  // combination: STRUCTURE <--> CONTACT
  else if ((qtype_1==NOX::NLN::StatusTest::quantity_structure and
        qtype_2==NOX::NLN::StatusTest::quantity_contact) or
       (qtype_1==NOX::NLN::StatusTest::quantity_contact and
        qtype_2==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_fres_contact_res_;
  // no combination was found
  else
    dserror("There is no combination type for the given quantity types! "
        "(quantity_1: %s, quantity_2: %s)",
        NOX::NLN::StatusTest::QuantityType2String(qtype_1).c_str(),
        NOX::NLN::StatusTest::QuantityType2String(qtype_2).c_str());

  return INPAR::STR::bop_and;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::BinaryOp STR::TIMINT::BaseDataSDyn::GetIncrComboType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  return GetIncrComboType(NOX::NLN::StatusTest::quantity_structure,qtype);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::BinaryOp STR::TIMINT::BaseDataSDyn::GetIncrComboType(
    const enum NOX::NLN::StatusTest::QuantityType& qtype_1,
    const enum NOX::NLN::StatusTest::QuantityType& qtype_2) const
{
  // combination: STRUCTURE <--> PRESSURE
  if ((qtype_1==NOX::NLN::StatusTest::quantity_structure and
       qtype_2==NOX::NLN::StatusTest::quantity_pressure) or
      (qtype_1==NOX::NLN::StatusTest::quantity_pressure and
       qtype_2==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_disp_pres_;
  // combination: STRUCTURE <--> EAS
  else if ((qtype_1==NOX::NLN::StatusTest::quantity_structure and
        qtype_2==NOX::NLN::StatusTest::quantity_eas) or
       (qtype_1==NOX::NLN::StatusTest::quantity_eas and
        qtype_2==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_disp_eas_incr_;
  // combination: STRUCTURE <--> PLASTICITY
  else if ((qtype_1==NOX::NLN::StatusTest::quantity_structure and
        qtype_2==NOX::NLN::StatusTest::quantity_plasticity) or
       (qtype_1==NOX::NLN::StatusTest::quantity_plasticity and
        qtype_2==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_disp_plast_incr_;
  // combination: STRUCTURE <--> CONTACT
  else if ((qtype_1==NOX::NLN::StatusTest::quantity_structure and
        qtype_2==NOX::NLN::StatusTest::quantity_contact) or
       (qtype_1==NOX::NLN::StatusTest::quantity_contact and
        qtype_2==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_disp_contact_lm_incr_;
  // no combination was found
  else
    dserror("There is no combination type for the given quantity types! "
        "(quantity_1: %s, quantity_2: %s)",
        NOX::NLN::StatusTest::QuantityType2String(qtype_1).c_str(),
        NOX::NLN::StatusTest::QuantityType2String(qtype_2).c_str());

  return INPAR::STR::bop_and;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::BinaryOp STR::TIMINT::BaseDataSDyn::GetResIncrComboType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype_res,
        const enum NOX::NLN::StatusTest::QuantityType& qtype_incr) const
{
  // combination: STRUCTURE (force/res) <--> STRUCTURE (displ/incr)
  if ((qtype_res==NOX::NLN::StatusTest::quantity_structure and
       qtype_incr==NOX::NLN::StatusTest::quantity_structure))
    return normcombo_fres_disp_;
  // no combination was found
  else
    dserror("There is no res-incr-combination type for the given quantity types! "
        "(quantity_res: %s, quantity_incr: %s)",
        NOX::NLN::StatusTest::QuantityType2String(qtype_res).c_str(),
        NOX::NLN::StatusTest::QuantityType2String(qtype_incr).c_str());

  return INPAR::STR::bop_and;
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
  dis_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  // velocities V_{n}
  vel_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  // accelerations A_{n}
  acc_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));

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
  modelevaluators_ = STR::MODELEVALUATOR::BuildModelEvaluators(DataSDyn().GetModelTypes());

  std::map<enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >::iterator me_iter;
  for (me_iter=ModelEvaluatorsPtr()->begin();me_iter!=ModelEvaluatorsPtr()->end();++me_iter)
  {
    me_iter->second->Init();
    me_iter->second->Setup();
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
  // we have to remove the const state, because the IO-routine expects
  // a non-const version.
  Teuchos::RCP<DRT::Discretization> actdis =
      Teuchos::rcp_const_cast<DRT::Discretization>(DataGlobalState().GetDiscret());
  IO::DiscretizationReader reader(actdis, step);
  if (step != reader.ReadInt("step"))
    dserror("Time step on file not equal to given step");

//  step_ = step;
//  stepn_ = step_ + 1;
//  time_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
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
