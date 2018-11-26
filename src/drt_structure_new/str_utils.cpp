/*-----------------------------------------------------------*/
/*!
\file str_utils.cpp

\brief Utility methods for the structural time integration.

\maintainer Matthias Mayr

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_integrator.H"
#include "str_model_evaluator_contact.H"
#include "str_model_evaluator_meshtying.H"
#include "str_model_evaluator_lagpenconstraint.H"
#include "str_nln_linearsystem_scaling.H"
#include "str_utils.H"

#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/lagpenconstraint_noxinterface.H"

#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_noxinterface.H"
#include "../drt_contact/meshtying_abstract_strategy.H"
#include "../drt_contact/meshtying_noxinterface.H"

#include "../drt_inpar/inpar_structure.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_structure_new/str_timint_basedatasdyn.H"

#include "../solver_nonlin_nox/nox_nln_aux.H"
#include "../solver_nonlin_nox/nox_nln_constraint_interface_preconditioner.H"
#include "../solver_nonlin_nox/nox_nln_constraint_interface_required.H"

#include "../linalg/linalg_utils.H"

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::ConditionNumber STR::NLN::Convert2NoxConditionNumberType(
    const INPAR::STR::ConditionNumber stype)
{
  switch (stype)
  {
    case INPAR::STR::ConditionNumber::max_min_ev_ratio:
      return NOX::NLN::LinSystem::ConditionNumber::max_min_ev_ratio;
    case INPAR::STR::ConditionNumber::one_norm:
      return NOX::NLN::LinSystem::ConditionNumber::one_norm;
    case INPAR::STR::ConditionNumber::inf_norm:
      return NOX::NLN::LinSystem::ConditionNumber::inf_norm;
    default:
      dserror("No known conversion.");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::Abstract::Vector::NormType STR::NLN::Convert2NoxNormType(
    const enum INPAR::STR::VectorNorm& normtype)
{
  enum NOX::Abstract::Vector::NormType nox_normtype = NOX::Abstract::Vector::TwoNorm;

  switch (normtype)
  {
    case INPAR::STR::norm_l2:
      nox_normtype = NOX::Abstract::Vector::TwoNorm;
      break;
    case INPAR::STR::norm_l1:
      nox_normtype = NOX::Abstract::Vector::OneNorm;
      break;
    case INPAR::STR::norm_inf:
      nox_normtype = NOX::Abstract::Vector::MaxNorm;
      break;
    case INPAR::STR::norm_rms:
    case INPAR::STR::norm_vague:
    default:
      dserror("Unknown conversion for the given vector norm type: \" %s \"!",
          INPAR::STR::VectorNormString(normtype).c_str());
      break;
  }  // switch case normtype

  return nox_normtype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::ConvertModelType2SolType(std::vector<enum NOX::NLN::SolutionType>& soltypes,
    std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& slinsolvers,
    const std::set<enum INPAR::STR::ModelType>& modeltypes,
    const std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver>>& mlinsolvers)
{
  // initialize the vector and/or force the length to zero
  if (soltypes.size() > 0)
  {
    soltypes.clear();
    slinsolvers.clear();
  }

  // pre-set the vector size
  soltypes.reserve(modeltypes.size());

  // The strings of the different enums have to fit!
  std::set<enum INPAR::STR::ModelType>::const_iterator mt_iter;
  for (mt_iter = modeltypes.begin(); mt_iter != modeltypes.end(); ++mt_iter)
  {
    const enum NOX::NLN::SolutionType soltype = ConvertModelType2SolType(*mt_iter);

    soltypes.push_back(soltype);
    // copy the linsolver pointers into the new map
    if (mlinsolvers.find(*mt_iter) != mlinsolvers.end())
      slinsolvers[soltype] = mlinsolvers.at(*mt_iter);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::SolutionType STR::NLN::ConvertModelType2SolType(
    const enum INPAR::STR::ModelType& modeltype, const bool& do_check)
{
  enum NOX::NLN::SolutionType soltype = NOX::NLN::sol_unknown;
  switch (modeltype)
  {
    case INPAR::STR::model_structure:
    case INPAR::STR::model_springdashpot:
    case INPAR::STR::model_monolithic_coupling:
    case INPAR::STR::model_partitioned_coupling:
    case INPAR::STR::model_beam_interaction_old:
    case INPAR::STR::model_browniandyn:
    case INPAR::STR::model_beaminteraction:
      soltype = NOX::NLN::sol_structure;
      break;
    case INPAR::STR::model_contact:
      soltype = NOX::NLN::sol_contact;
      break;
    case INPAR::STR::model_meshtying:
      soltype = NOX::NLN::sol_meshtying;
      break;
    case INPAR::STR::model_cardiovascular0d:
      soltype = NOX::NLN::sol_cardiovascular0d;
      break;
    case INPAR::STR::model_lag_pen_constraint:
      soltype = NOX::NLN::sol_lag_pen_constraint;
      break;
    default:
      // check if the corresponding enum could be found.
      if (do_check)
        dserror(
            "The corresponding solution-type was not found. "
            "Given string: %s",
            INPAR::STR::ModelTypeString(modeltype).c_str());
      break;
  }

  return soltype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::ModelType STR::NLN::ConvertSolType2ModelType(
    const enum NOX::NLN::SolutionType& soltype, const bool& do_check)
{
  enum INPAR::STR::ModelType modeltype = INPAR::STR::model_vague;
  switch (soltype)
  {
    case NOX::NLN::sol_structure:
      modeltype = INPAR::STR::model_structure;
      break;
    case NOX::NLN::sol_contact:
      modeltype = INPAR::STR::model_contact;
      break;
    case NOX::NLN::sol_meshtying:
      modeltype = INPAR::STR::model_meshtying;
      break;
    case NOX::NLN::sol_cardiovascular0d:
      modeltype = INPAR::STR::model_cardiovascular0d;
      break;
    case NOX::NLN::sol_lag_pen_constraint:
      modeltype = INPAR::STR::model_lag_pen_constraint;
      break;
    default:
      // check if the corresponding enum could be found.
      if (do_check)
        dserror(
            "The corresponding model-type was not found. "
            "Given string: %s",
            NOX::NLN::SolutionType2String(soltype).c_str());
      break;
  }

  return modeltype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::ModelType STR::NLN::ConvertQuantityType2ModelType(
    const enum NOX::NLN::StatusTest::QuantityType& qtype, const bool& do_check)
{
  const NOX::NLN::SolutionType st = NOX::NLN::AUX::ConvertQuantityType2SolutionType(qtype);
  return ConvertSolType2ModelType(st, do_check);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::EleTech STR::NLN::ConvertQuantityType2EleTech(
    const enum NOX::NLN::StatusTest::QuantityType& qtype)
{
  enum INPAR::STR::EleTech eletech = INPAR::STR::eletech_pressure;
  switch (qtype)
  {
    case NOX::NLN::StatusTest::quantity_pressure:
      eletech = INPAR::STR::eletech_pressure;
      break;
    case NOX::NLN::StatusTest::quantity_plasticity:
      eletech = INPAR::STR::eletech_plasticity;
      break;
    case NOX::NLN::StatusTest::quantity_eas:
      eletech = INPAR::STR::eletech_eas;
      break;
    default:
      dserror("Cannot convert QuantityType %s to EleTech.",
          NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      break;
  }

  return eletech;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::OptimizationProblemType STR::NLN::OptimizationType(
    const std::vector<enum NOX::NLN::SolutionType>& soltypes)
{
  enum NOX::NLN::OptimizationProblemType opttype = NOX::NLN::opt_unconstrained;
  std::vector<enum NOX::NLN::SolutionType>::const_iterator st_iter;

  for (st_iter = soltypes.begin(); st_iter != soltypes.end(); ++st_iter)
  {
    switch (*st_iter)
    {
      // -----------------------------------
      // Inequality constraint
      // active set decision necessary
      // saddle point structure or condensed
      // -----------------------------------
      case NOX::NLN::sol_contact:
        return NOX::NLN::opt_inequality_constrained;
        break;
      // -----------------------------------
      // Equality constraint
      // no active set decision necessary
      // saddle point structure or condensed
      // -----------------------------------
      case NOX::NLN::sol_meshtying:
      case NOX::NLN::sol_lag_pen_constraint:
        opttype = NOX::NLN::opt_equality_constrained;
        break;
      // -----------------------------------
      // Unconstrained problem
      // pure structural problem
      // no saddle point structure
      // -----------------------------------
      case NOX::NLN::sol_structure:
      case NOX::NLN::sol_cardiovascular0d:
      default:
        break;
    }
  }

  return opttype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::CreateConstraintInterfaces(NOX::NLN::CONSTRAINT::ReqInterfaceMap& iconstr,
    STR::Integrator& integrator, const std::vector<enum NOX::NLN::SolutionType>& soltypes)
{
  if (iconstr.size() > 0) iconstr.clear();

  std::vector<enum NOX::NLN::SolutionType>::const_iterator st_iter;
  for (st_iter = soltypes.begin(); st_iter != soltypes.end(); ++st_iter)
  {
    switch (*st_iter)
    {
      case NOX::NLN::sol_contact:
      {
        STR::MODELEVALUATOR::Generic& model = integrator.Evaluator(INPAR::STR::model_contact);
        STR::MODELEVALUATOR::Contact& contact_model =
            dynamic_cast<STR::MODELEVALUATOR::Contact&>(model);
        iconstr[NOX::NLN::sol_contact] = contact_model.StrategyPtr()->NoxInterfacePtr();
        break;
      }
      case NOX::NLN::sol_meshtying:
      {
        STR::MODELEVALUATOR::Generic& model = integrator.Evaluator(INPAR::STR::model_meshtying);
        STR::MODELEVALUATOR::Meshtying& mt_model =
            dynamic_cast<STR::MODELEVALUATOR::Meshtying&>(model);
        iconstr[NOX::NLN::sol_meshtying] = mt_model.StrategyPtr()->NoxInterfacePtr();
        break;
      }
      case NOX::NLN::sol_lag_pen_constraint:
      {
        STR::MODELEVALUATOR::Generic& model =
            integrator.Evaluator(INPAR::STR::model_lag_pen_constraint);
        STR::MODELEVALUATOR::LagPenConstraint& lagpenconstraint_model =
            dynamic_cast<STR::MODELEVALUATOR::LagPenConstraint&>(model);
        iconstr[NOX::NLN::sol_lag_pen_constraint] = lagpenconstraint_model.NoxInterfacePtr();
        break;
      }
      default:
        break;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::CreateConstraintPreconditioner(NOX::NLN::CONSTRAINT::PrecInterfaceMap& iconstr_prec,
    STR::Integrator& integrator, const std::vector<enum NOX::NLN::SolutionType>& soltypes)
{
  if (iconstr_prec.size() > 0) iconstr_prec.clear();

  std::vector<enum NOX::NLN::SolutionType>::const_iterator st_iter;
  for (st_iter = soltypes.begin(); st_iter != soltypes.end(); ++st_iter)
  {
    switch (*st_iter)
    {
      case NOX::NLN::sol_contact:
      {
        STR::MODELEVALUATOR::Generic& model = integrator.Evaluator(INPAR::STR::model_contact);
        STR::MODELEVALUATOR::Contact& contact_model =
            dynamic_cast<STR::MODELEVALUATOR::Contact&>(model);
        /* Actually we use the underlying MORTAR::StrategyBase as Preconditioner
         * interface. Nevertheless, the implementations can differ for the
         * contact/meshtying cases. */
        iconstr_prec[NOX::NLN::sol_contact] = contact_model.StrategyPtr();
        break;
      }
      case NOX::NLN::sol_meshtying:
      {
        STR::MODELEVALUATOR::Generic& model = integrator.Evaluator(INPAR::STR::model_meshtying);
        STR::MODELEVALUATOR::Meshtying& mt_model =
            dynamic_cast<STR::MODELEVALUATOR::Meshtying&>(model);
        iconstr_prec[NOX::NLN::sol_meshtying] = mt_model.StrategyPtr();
        break;
      }
      case NOX::NLN::sol_lag_pen_constraint:
      {
        STR::MODELEVALUATOR::Generic& model =
            integrator.Evaluator(INPAR::STR::model_lag_pen_constraint);
        STR::MODELEVALUATOR::LagPenConstraint& lagpenconstraint_model =
            dynamic_cast<STR::MODELEVALUATOR::LagPenConstraint&>(model);
        iconstr_prec[NOX::NLN::sol_lag_pen_constraint] =
            lagpenconstraint_model.NoxInterfacePrecPtr();
        break;
      }
      default:
        // do nothing
        break;
    }  // switch (*st_iter)
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::CreateScaling(Teuchos::RCP<NOX::Epetra::Scaling>& iscale,
    STR::TIMINT::BaseDataSDyn DataSDyn, STR::TIMINT::BaseDataGlobalState GState)
{
  if (DataSDyn.GetSTCAlgoType() != INPAR::STR::stc_none)
    iscale = Teuchos::rcp(new STR::NLN::LinSystem::StcScaling(DataSDyn, GState));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::GetTimIntFactor()
{
  const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
  const INPAR::STR::DynamicType dyntype =
      DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP");

  double timintfactor = 0.0;

  switch (dyntype)
  {
    case INPAR::STR::dyna_statics:
    {
      timintfactor = 0.0;
      break;
    }
    case INPAR::STR::dyna_genalpha:
    {
      timintfactor = sdynparams.sublist("GENALPHA").get<double>("ALPHA_F");
      break;
    }
    case INPAR::STR::dyna_gemm:
    {
      timintfactor = sdynparams.sublist("GEMM").get<double>("ALPHA_F");
      break;
    }
    case INPAR::STR::dyna_onesteptheta:
    {
      timintfactor = 1.0 - sdynparams.sublist("ONESTEPTHETA").get<double>("THETA");
      break;
    }
    default:
    {
      dserror("Time integration factor not been set for time integration scheme '%s'",
          INPAR::STR::DynamicTypeString(dyntype));
      break;
    }
  }

  return timintfactor;
}
