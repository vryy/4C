/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a mixture growth strategy interface

\level 3
*/
/*----------------------------------------------------------------------*/
#include "mixture_growth_strategy.H"
#include "globalproblem.H"
#include "matpar_bundle.H"
#include "material_service.H"
#include "mixture_growth_strategy_anisotropic.H"
#include "mixture_growth_strategy_isotropic.H"
#include "mixture_growth_strategy_stiffness.H"

MIXTURE::PAR::MixtureGrowthStrategy::MixtureGrowthStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MAT::PAR::Parameter(matdata)
{
}

MIXTURE::PAR::MixtureGrowthStrategy* MIXTURE::PAR::MixtureGrowthStrategy::Factory(const int matid)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
  {
    dserror("List of materials cannot be accessed in the global problem instance.");
  }

  // yet another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
  {
    dserror("List of materials in the global problem instance is empty.");
  }

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matid);

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
      dserror("The referenced material with id %d is not registered as a mixture growth strategy!",
          matid);
  }
  return nullptr;
}