/*----------------------------------------------------------------------*/
/*! \file

\brief General prestress strategy for mixture constituents

\level 3


*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_prestress_strategy.hpp"

#include "baci_global_data.hpp"
#include "baci_inpar_material.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_par_material.hpp"
#include "baci_mat_service.hpp"
#include "baci_mixture_prestress_strategy_constant.hpp"
#include "baci_mixture_prestress_strategy_isocyl.hpp"
#include "baci_mixture_prestress_strategy_iterative.hpp"
#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

// Prestress stragegy factory generates the prestress strategy for a specific material id
MIXTURE::PAR::PrestressStrategy* MIXTURE::PAR::PrestressStrategy::Factory(int matid)
{
  // for the sake of safety
  if (GLOBAL::Problem::Instance()->Materials() == Teuchos::null)
  {
    dserror("List of materials cannot be accessed in the global problem instance.");
  }

  // yet another safety check
  if (GLOBAL::Problem::Instance()->Materials()->Num() == 0)
  {
    dserror("List of materials in the global problem instance is empty.");
  }

  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(matid);

  switch (curmat->Type())
  {
    case INPAR::MAT::mix_prestress_strategy_cylinder:
    {
      return MAT::CreateMaterialParameterInstance<MIXTURE::PAR::IsotropicCylinderPrestressStrategy>(
          curmat);
    }
    case INPAR::MAT::mix_prestress_strategy_iterative:
    {
      return MAT::CreateMaterialParameterInstance<MIXTURE::PAR::IterativePrestressStrategy>(curmat);
    }
    case INPAR::MAT::mix_prestress_strategy_constant:
    {
      return MAT::CreateMaterialParameterInstance<MIXTURE::PAR::ConstantPrestressStrategy>(curmat);
    }
    default:
      dserror(
          "The referenced material with id %d is not registered as a prestress strategy!", matid);
  }

  return nullptr;
}
BACI_NAMESPACE_CLOSE
