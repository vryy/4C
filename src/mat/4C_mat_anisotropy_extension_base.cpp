/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a base anisotropy extension to be used by anisotropic materials with
@Mat::Anisotropy

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_mat_anisotropy_extension_base.hpp"

#include "4C_mat_anisotropy.hpp"

FOUR_C_NAMESPACE_OPEN

void Mat::BaseAnisotropyExtension::set_anisotropy(Mat::Anisotropy& anisotropy)
{
  anisotropy_ = Teuchos::rcpFromRef(anisotropy);
}

FOUR_C_NAMESPACE_CLOSE
