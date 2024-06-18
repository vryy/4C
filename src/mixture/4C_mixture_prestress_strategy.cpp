/*----------------------------------------------------------------------*/
/*! \file

\brief General prestress strategy for mixture constituents

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_prestress_strategy.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_material.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_prestress_strategy_constant.hpp"
#include "4C_mixture_prestress_strategy_isocyl.hpp"
#include "4C_mixture_prestress_strategy_iterative.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

// Prestress stragegy factory generates the prestress strategy for a specific material id
MIXTURE::PAR::PrestressStrategy* MIXTURE::PAR::PrestressStrategy::factory(int matid)
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
    case Core::Materials::mix_prestress_strategy_cylinder:
    {
      return Mat::create_material_parameter_instance<
          MIXTURE::PAR::IsotropicCylinderPrestressStrategy>(curmat);
    }
    case Core::Materials::mix_prestress_strategy_iterative:
    {
      return Mat::create_material_parameter_instance<MIXTURE::PAR::IterativePrestressStrategy>(
          curmat);
    }
    case Core::Materials::mix_prestress_strategy_constant:
    {
      return Mat::create_material_parameter_instance<MIXTURE::PAR::ConstantPrestressStrategy>(
          curmat);
    }
    default:
      FOUR_C_THROW(
          "The referenced material with id %d is not registered as a prestress strategy!", matid);
  }

  return nullptr;
}
FOUR_C_NAMESPACE_CLOSE
