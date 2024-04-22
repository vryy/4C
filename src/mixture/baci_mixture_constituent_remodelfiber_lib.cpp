/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of helper functions for the remodel fiber constituent
\level 3
*/
/*----------------------------------------------------------------------*/

#include "baci_mixture_constituent_remodelfiber_lib.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_service.hpp"
#include "baci_mixture_constituent_remodelfiber_material.hpp"
#include "baci_mixture_constituent_remodelfiber_material_exponential.hpp"
#include "baci_mixture_constituent_remodelfiber_material_exponential_active.hpp"

FOUR_C_NAMESPACE_OPEN

[[nodiscard]] const MIXTURE::PAR::RemodelFiberMaterial<double>* MIXTURE::PAR::FiberMaterialFactory(
    int matid)
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
    case INPAR::MAT::mix_remodelfiber_material_exponential:
      return MAT::CreateMaterialParameterInstance<
          MIXTURE::PAR::RemodelFiberMaterialExponential<double>>(curmat);
    case INPAR::MAT::mix_remodelfiber_material_exponential_active:
      return MAT::CreateMaterialParameterInstance<
          MIXTURE::PAR::RemodelFiberMaterialExponentialActive<double>>(curmat);
    default:
      FOUR_C_THROW(
          "The referenced material with id %d is not registered as a remodel fiber material!",
          matid);
  }

  // we will not end up here, so make the compiler happy
  std::exit(1);
}
FOUR_C_NAMESPACE_CLOSE
