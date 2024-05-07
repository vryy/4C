/*----------------------------------------------------------------------*/
/*! \file
\brief A material parameter container

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_material_parameter_base.hpp"

#include "4C_material_input_base.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


CORE::MAT::PAR::Parameter::Parameter(Teuchos::RCP<const CORE::MAT::PAR::Material> matdata)
    : id_(matdata->Id()), type_(matdata->Type()), name_(matdata->Name())
{
}

void CORE::MAT::PAR::Parameter::SetParameter(
    int parametername, Teuchos::RCP<Epetra_Vector> myparameter)
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

void CORE::MAT::PAR::Parameter::SetParameter(int parametername, const double val, const int eleGID)
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

double CORE::MAT::PAR::Parameter::GetParameter(int parametername, const int EleId)
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


FOUR_C_NAMESPACE_CLOSE
