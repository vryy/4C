/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for plasticity
\level 2
*/

/*----------------------------------------------------------------------*/


#include "4C_inpar_plasticity.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_inpar_tsi.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::Plasticity::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*----------------------------------------------------------------------*/
  /* parameters for semi-smooth Newton plasticity algorithm */
  Teuchos::ParameterList& iplast = list->sublist("SEMI-SMOOTH PLASTICITY", false, "");

  Core::UTILS::DoubleParameter(
      "SEMI_SMOOTH_CPL", 1.0, "Weighting factor cpl for semi-smooth PDASS", &iplast);
  Core::UTILS::DoubleParameter(
      "STABILIZATION_S", 1.0, "Stabilization factor s for semi-smooth PDASS", &iplast);

  // solver convergence test parameters for semi-smooth plasticity formulation
  setStringToIntegralParameter<int>("NORMCOMBI_RESFPLASTCONSTR", "And",
      "binary operator to combine plasticity constraints and residual force values",
      tuple<std::string>("And", "Or"), tuple<int>(Inpar::STR::bop_and, Inpar::STR::bop_or),
      &iplast);

  setStringToIntegralParameter<int>("NORMCOMBI_DISPPLASTINCR", "And",
      "binary operator to combine displacement increments and plastic flow (Delta Lp) increment "
      "values",
      tuple<std::string>("And", "Or"), tuple<int>(Inpar::STR::bop_and, Inpar::STR::bop_or),
      &iplast);

  Core::UTILS::DoubleParameter("TOLPLASTCONSTR", 1.0E-8,
      "tolerance in the plastic constraint norm for the newton iteration", &iplast);
  Core::UTILS::DoubleParameter("TOLDELTALP", 1.0E-8,
      "tolerance in the plastic flow (Delta Lp) norm for the Newton iteration", &iplast);

  setStringToIntegralParameter<int>("NORMCOMBI_EASRES", "And",
      "binary operator to combine EAS-residual and residual force values",
      tuple<std::string>("And", "Or"), tuple<int>(Inpar::STR::bop_and, Inpar::STR::bop_or),
      &iplast);

  setStringToIntegralParameter<int>("NORMCOMBI_EASINCR", "And",
      "binary operator to combine displacement increments and EAS increment values",
      tuple<std::string>("And", "Or"), tuple<int>(Inpar::STR::bop_and, Inpar::STR::bop_or),
      &iplast);

  Core::UTILS::DoubleParameter(
      "TOLEASRES", 1.0E-8, "tolerance in the EAS residual norm for the newton iteration", &iplast);
  Core::UTILS::DoubleParameter("TOLEASINCR", 1.0E-8,
      "tolerance in the EAS increment norm for the Newton iteration", &iplast);

  setStringToIntegralParameter<int>("DISSIPATION_MODE", "pl_multiplier",
      "method to calculate the plastic dissipation",
      tuple<std::string>("pl_multiplier", "pl_flow", "Taylor_Quinney"),
      tuple<int>(Inpar::TSI::pl_multiplier, Inpar::TSI::pl_flow, Inpar::TSI::Taylor_Quinney),
      &iplast);
}

FOUR_C_NAMESPACE_CLOSE
