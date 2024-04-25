/*----------------------------------------------------------------------*/
/*! \file
\brief A material parameter container

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_par_parameter.hpp"

#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_par_material.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


MAT::PAR::Parameter::Parameter(Teuchos::RCP<const MAT::PAR::Material> matdata)
    : id_(matdata->Id()), type_(matdata->Type()), name_(matdata->Name())
{
}

void MAT::PAR::Parameter::SetParameter(int parametername, Teuchos::RCP<Epetra_Vector> myparameter)
{
  // security check
  int length_old = matparams_.at(parametername)->GlobalLength();
  int length_new = myparameter->GlobalLength();

  // we had spatially constant matparameter and want elementwise, thats ok
  if (length_old == 1 && (length_old < length_new))
  {
    matparams_.at(parametername) = Teuchos::rcp(new Epetra_Vector(*myparameter));
  }
  // we had spatially constant matparameter and want constant, thats ok
  else if (length_old == 1 && (length_old == length_new))
  {
    matparams_.at(parametername) = Teuchos::rcp(new Epetra_Vector(*myparameter));
  }
  // we had elementwise matparameter and want elementwise, thats ok
  else if (length_old > 1 && (length_old == length_new))
  {
    matparams_.at(parametername) = Teuchos::rcp(new Epetra_Vector(*myparameter));
  }
  // we had elementwise matparameter and want constant, thats OK, simply set same value in all
  // entries
  else if (length_old > 1 && (length_new == 1))
  {
    matparams_.at(parametername)->PutScalar((*myparameter)[0]);
  }
  // we had elementwise matparameter and elementwise but there is  size mismatch, thats  not ok
  else if (length_old > 1 && (length_old != length_new))
  {
    FOUR_C_THROW("Problem setting elementwise material parameter: Size mismatch ");
  }
  else
  {
    FOUR_C_THROW("Can not set material parameter: Unknown Problem");
  }
}

void MAT::PAR::Parameter::SetParameter(int parametername, const double val, const int eleGID)
{
  // check if we own this element
  Teuchos::RCP<Epetra_Vector> fool = matparams_.at(parametername);
  if (matparams_.at(parametername)->GlobalLength() == 1)
    FOUR_C_THROW("spatially constant parameters! Cannot set elementwise parameter!");

  if (!matparams_.at(parametername)->Map().MyGID(eleGID))
    FOUR_C_THROW("I do not have this element");

  // otherwise set parameter for element
  // calculate LID here, instead of before each call
  (*matparams_.at(parametername))[matparams_[parametername]->Map().LID(eleGID)] = val;
}

void MAT::PAR::Parameter::ExpandParametersToEleColLayout()
{
  for (auto& matparam : matparams_)
  {
    // only do this for vectors with one entry
    if (matparam->GlobalLength() == 1)
    {
      // get value of element
      double temp = (*matparam)[0];
      // put new RCP<Epetra_Vector> in matparams struct
      Teuchos::RCP<Epetra_Vector> temp2 = Teuchos::rcp(new Epetra_Vector(
          *(GLOBAL::Problem::Instance()->GetDis("structure")->ElementColMap()), true));
      temp2->PutScalar(temp);
      matparam = temp2;
    }
  }
}
double MAT::PAR::Parameter::GetParameter(int parametername, const int EleId)
{
  // check if we have an element based value via size
  if (matparams_[parametername]->GlobalLength() == 1)
  {
    // we have a global value hence we directly return the first entry
    return (*matparams_[parametername])[0];
  }
  // If someone calls this functions without a valid EleID and we have element based values throw
  // error
  else if (EleId < 0 && matparams_[parametername]->GlobalLength() > 1)
  {
    FOUR_C_THROW("Global mat parameter requested but we have elementwise mat params");
    return 0.0;
  }
  // otherwise just return the element specific value
  else
  {
    // calculate LID here, instead of before each call
    return (*matparams_[parametername])[matparams_[parametername]->Map().LID(EleId)];
  }
}

MAT::PAR::ParameterAniso::ParameterAniso(Teuchos::RCP<const MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  // get MAT ID for definiton of structural tensor
  int mat_id_structural_tensor = *matdata->Get<int>("STR_TENS_ID");
  // get pointer to material
  Teuchos::RCP<MAT::PAR::Material> mat_str_tens =
      GLOBAL::Problem::Instance()->Materials()->ById(mat_id_structural_tensor);
  // construct parameter class
  if (mat_str_tens->Parameter() == nullptr)
    mat_str_tens->SetParameter(new MAT::ELASTIC::PAR::StructuralTensorParameter(mat_str_tens));
  auto* params =
      static_cast<MAT::ELASTIC::PAR::StructuralTensorParameter*>(mat_str_tens->Parameter());
  // get type of strategy
  std::string strategy = *mat_str_tens->Get<std::string>("STRATEGY");

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
