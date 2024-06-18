/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a mixture growth strategy interface

\level 3
*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_growth_strategy.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_mixture_growth_strategy_anisotropic.hpp"
#include "4C_mixture_growth_strategy_isotropic.hpp"
#include "4C_mixture_growth_strategy_stiffness.hpp"

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::MixtureGrowthStrategy::MixtureGrowthStrategy(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Core::Mat::PAR::Parameter(matdata)
{
}

MIXTURE::PAR::MixtureGrowthStrategy* MIXTURE::PAR::MixtureGrowthStrategy::Factory(const int matid)
{
  // for the sake of safety
  if (Global::Problem::Instance()->Materials() == Teuchos::null)
  {
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  }

  // yet another safety check
  if (Global::Problem::Instance()->Materials()->Num() == 0)
  {
    FOUR_C_THROW("List of materials in the global problem instance is empty.");
  }

  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);

  switch (curmat->Type())
  {
    case Core::Materials::mix_growth_strategy_isotropic:
    {
      return Mat::create_material_parameter_instance<MIXTURE::PAR::IsotropicGrowthStrategy>(curmat);
    }
    case Core::Materials::mix_growth_strategy_anisotropic:
    {
      return Mat::create_material_parameter_instance<MIXTURE::PAR::AnisotropicGrowthStrategy>(
          curmat);
    }
    case Core::Materials::mix_growth_strategy_stiffness:
    {
      return Mat::create_material_parameter_instance<MIXTURE::PAR::StiffnessGrowthStrategy>(curmat);
    }
    default:
      FOUR_C_THROW(
          "The referenced material with id %d is not registered as a mixture growth strategy!",
          matid);
  }
  return nullptr;
}
FOUR_C_NAMESPACE_CLOSE
