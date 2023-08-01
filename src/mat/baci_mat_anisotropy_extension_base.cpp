/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a base anisotropy extension to be used by anisotropic materials with
@MAT::Anisotropy

\level 3


*/
/*----------------------------------------------------------------------*/

#include "baci_mat_anisotropy_extension_base.H"

#include "baci_mat_anisotropy.H"

void MAT::BaseAnisotropyExtension::SetAnisotropy(MAT::Anisotropy& anisotropy)
{
  anisotropy_ = Teuchos::rcpFromRef(anisotropy);
}
