/*----------------------------------------------------------------------*/
/*! \file
\brief A material parameter container

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_material_parameter_base.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


Core::Mat::PAR::Parameter::Parameter(Data data) : data_(std::move(data)) {}

double Core::Mat::PAR::Parameter::GetParameter(int parametername, const int EleId)
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
