/*----------------------------------------------------------------------*/
/*!
\brief
Former file of Lena Yoshihara

\level 1

<pre>
\maintainer Fabian Braeu
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "matpar_parameter.H"
#include "matpar_material.H"
#include "matpar_bundle.H"

#include "../drt_matelast/elast_aniso_structuraltensor_strategy.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io_pstream.H"
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
    dserror("Problem setting elementwise material parameter: Size mismatch ");
  }
  else
  {
    dserror("Can not set material parameter: Unknown Problem");
  }
}

void MAT::PAR::Parameter::SetParameter(int parametername, const double val, const int eleGID)
{
  // check if we own this element
  Teuchos::RCP<Epetra_Vector> fool = matparams_.at(parametername);
  if (matparams_.at(parametername)->GlobalLength() == 1)
    dserror("spatially constant parameters! Cannot set elementwise parameter!");

  if (!matparams_.at(parametername)->Map().MyGID(eleGID)) dserror("I do not have this element");

  // otherwise set parameter for element
  // calculate LID here, instead of before each call              01/2017 birzle
  (*matparams_.at(parametername))[matparams_[parametername]->Map().LID(eleGID)] = val;
  // old: (*matparams_.at(parametername))[eleLID]=val;
}

void MAT::PAR::Parameter::ExpandParametersToEleColLayout()
{
  for (unsigned int i = 0; i < matparams_.size(); i++)
  {
    // only do this for vectors with one entry
    if (matparams_.at(i)->GlobalLength() == 1)
    {
      // get value of element
      double temp = (*matparams_.at(i))[0];
      // put new RCP<Epetra_Vector> in matparams struct
      Teuchos::RCP<Epetra_Vector> temp2 = Teuchos::rcp(new Epetra_Vector(
          *(DRT::Problem::Instance()->GetDis("structure")->ElementColMap()), true));
      temp2->PutScalar(temp);
      matparams_.at(i) = temp2;
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
    dserror("Global mat parameter requested but we have elementwise mat params");
    return 0.0;
  }
  // otherwise just return the element specific value
  else
  {
    // calculate LID here, instead of before each call              01/2017 birzle
    return (*matparams_[parametername])[matparams_[parametername]->Map().LID(EleId)];
  }
}
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ParameterAniso::ParameterAniso(Teuchos::RCP<const MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  // get MAT ID for definiton of structural tensor
  int mat_id_structural_tensor = matdata->GetInt("STR_TENS_ID");
  // get pointer to material
  Teuchos::RCP<MAT::PAR::Material> mat_str_tens =
      DRT::Problem::Instance()->Materials()->ById(mat_id_structural_tensor);
  // construct parameter class
  if (mat_str_tens->Parameter() == NULL)
    mat_str_tens->SetParameter(new MAT::ELASTIC::PAR::StructuralTensorParameter(mat_str_tens));
  MAT::ELASTIC::PAR::StructuralTensorParameter* params =
      static_cast<MAT::ELASTIC::PAR::StructuralTensorParameter*>(mat_str_tens->Parameter());
  // get type of strategy
  std::string strategy = *mat_str_tens->Get<std::string>("STRATEGY");

  // construct strategy
  if (strategy == "Standard")
    structural_tensor_strategy_ =
        Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(params));
  else if (strategy == "ByDistributionFunction")
    structural_tensor_strategy_ =
        Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyByDistributionFunction(params));
  else if (strategy == "DispersedTransverselyIsotropic")
    structural_tensor_strategy_ = Teuchos::rcp(
        new MAT::ELASTIC::StructuralTensorStrategyDispersedTransverselyIsotropic(params));
  else
    dserror("Unknown type of structural tensor strategy for anisotropic material chosen.");
}
/*----------------------------------------------------------------------*/
