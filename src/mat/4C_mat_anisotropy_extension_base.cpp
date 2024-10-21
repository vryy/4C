#include "4C_mat_anisotropy_extension_base.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_mat_anisotropy.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

void Mat::BaseAnisotropyExtension::set_anisotropy(Mat::Anisotropy& anisotropy)
{
  anisotropy_ = Teuchos::rcpFromRef(anisotropy);
}

FOUR_C_NAMESPACE_CLOSE
