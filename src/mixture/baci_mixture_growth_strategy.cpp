/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a mixture growth strategy interface

\level 3
*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_growth_strategy.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_service.hpp"
#include "baci_mixture_growth_strategy_anisotropic.hpp"
#include "baci_mixture_growth_strategy_isotropic.hpp"
#include "baci_mixture_growth_strategy_stiffness.hpp"

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::MixtureGrowthStrategy::MixtureGrowthStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MAT::PAR::Parameter(matdata)
{
}

MIXTURE::PAR::MixtureGrowthStrategy* MIXTURE::PAR::MixtureGrowthStrategy::Factory(const int matid)
{
  // for the sake of safety
  if (GLOBAL::Problem::Instance()->Materials() == Teuchos::null)
  {
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  }

  // yet another safety check
  if (GLOBAL::Problem::Instance()->Materials()->Num() == 0)
  {
    FOUR_C_THROW("List of materials in the global problem instance is empty.");
  }

  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(matid);

  switch (curmat->Type())
  {
    case INPAR::MAT::mix_growth_strategy_isotropic:
    {
      return MAT::CreateMaterialParameterInstance<MIXTURE::PAR::IsotropicGrowthStrategy>(curmat);
    }
    case INPAR::MAT::mix_growth_strategy_anisotropic:
    {
      return MAT::CreateMaterialParameterInstance<MIXTURE::PAR::AnisotropicGrowthStrategy>(curmat);
    }
    case INPAR::MAT::mix_growth_strategy_stiffness:
    {
      return MAT::CreateMaterialParameterInstance<MIXTURE::PAR::StiffnessGrowthStrategy>(curmat);
    }
    default:
      FOUR_C_THROW(
          "The referenced material with id %d is not registered as a mixture growth strategy!",
          matid);
  }
  return nullptr;
}
FOUR_C_NAMESPACE_CLOSE
