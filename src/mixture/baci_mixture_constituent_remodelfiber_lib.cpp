/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of helper functions for the remodel fiber constituent
\level 3
*/
/*----------------------------------------------------------------------*/

#include "baci_mixture_constituent_remodelfiber_lib.H"

#include "baci_global_data.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_service.H"
#include "baci_mixture_constituent_remodelfiber_material.H"
#include "baci_mixture_constituent_remodelfiber_material_exponential.H"
#include "baci_mixture_constituent_remodelfiber_material_exponential_active.H"

BACI_NAMESPACE_OPEN

[[nodiscard]] const MIXTURE::PAR::RemodelFiberMaterial<double>* MIXTURE::PAR::FiberMaterialFactory(
    int matid)
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
    case INPAR::MAT::mix_remodelfiber_material_exponential:
      return MAT::CreateMaterialParameterInstance<
          MIXTURE::PAR::RemodelFiberMaterialExponential<double>>(curmat);
    case INPAR::MAT::mix_remodelfiber_material_exponential_active:
      return MAT::CreateMaterialParameterInstance<
          MIXTURE::PAR::RemodelFiberMaterialExponentialActive<double>>(curmat);
    default:
      dserror("The referenced material with id %d is not registered as a remodel fiber material!",
          matid);
  }

  // we will not end up here, so make the compiler happy
  std::exit(1);
}
BACI_NAMESPACE_CLOSE
