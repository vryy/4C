/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a base class for anisotropic material parameters

\level 1

*/

/*----------------------------------------------------------------------*/
/* macros */


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_mat_par_aniso.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"


FOUR_C_NAMESPACE_OPEN


MAT::PAR::ParameterAniso::ParameterAniso(Teuchos::RCP<const CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  // get MAT ID for definiton of structural tensor
  int mat_id_structural_tensor = matdata->Get<int>("STR_TENS_ID");
  // get pointer to material
  Teuchos::RCP<CORE::MAT::PAR::Material> mat_str_tens =
      GLOBAL::Problem::Instance()->Materials()->ById(mat_id_structural_tensor);
  // construct parameter class
  if (mat_str_tens->Parameter() == nullptr)
    mat_str_tens->SetParameter(new MAT::ELASTIC::PAR::StructuralTensorParameter(mat_str_tens));
  auto* params =
      static_cast<MAT::ELASTIC::PAR::StructuralTensorParameter*>(mat_str_tens->Parameter());
  // get type of strategy
  std::string strategy = mat_str_tens->Get<std::string>("STRATEGY");

  // construct strategy
  if (strategy == "Standard")
  {
    structural_tensor_strategy_ =
        Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(params));
  }
  else if (strategy == "ByDistributionFunction")
  {
    structural_tensor_strategy_ =
        Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyByDistributionFunction(params));
  }
  else if (strategy == "DispersedTransverselyIsotropic")
  {
    structural_tensor_strategy_ = Teuchos::rcp(
        new MAT::ELASTIC::StructuralTensorStrategyDispersedTransverselyIsotropic(params));
  }
  else
    FOUR_C_THROW("Unknown type of structural tensor strategy for anisotropic material chosen.");
}

FOUR_C_NAMESPACE_CLOSE
